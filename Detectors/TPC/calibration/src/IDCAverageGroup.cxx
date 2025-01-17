// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TPCCalibration/IDCAverageGroup.h"
#include "TPCCalibration/IDCGroup.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "TPCCalibration/IDCGroupingParameter.h"
#include "TPCBase/Mapper.h"

#include "TFile.h"
#include "TKey.h"
#include "TPCBase/Painter.h"
#include "TH2Poly.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TKey.h"
#include "Framework/Logger.h"

#if (defined(WITH_OPENMP) || defined(_OPENMP)) && !defined(__CLING__)
#include <omp.h>
#else
static inline int omp_get_thread_num() { return 0; }
#endif

o2::tpc::IDCAverageGroup::IDCAverageGroup(const unsigned char groupPads, const unsigned char groupRows, const unsigned char groupLastRowsThreshold, const unsigned char groupLastPadsThreshold, const unsigned int region, const Sector sector, const float sigma)
  : mIDCsGrouped{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, region}, mSector{sector}, mSigma{sigma}, mRobustAverage(sNThreads)
{
  unsigned int maxValues = 0;
  for (unsigned int i = 0; i < Mapper::NREGIONS; ++i) {
    const unsigned int maxGroup = (mIDCsGrouped.getGroupRows() + mIDCsGrouped.getGroupLastRowsThreshold()) * (mIDCsGrouped.getGroupPads() + mIDCsGrouped.getGroupLastPadsThreshold() + Mapper::ADDITIONALPADSPERROW[i].back());
    if (maxGroup > maxValues) {
      maxValues = maxGroup;
    }
  }

  for (auto& rob : mRobustAverage) {
    rob.reserve(maxValues);
  }
}

void o2::tpc::IDCAverageGroup::processIDCs()
{
  const static auto& paramIDCGroup = ParameterIDCGroup::Instance();

#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int integrationInterval = 0; integrationInterval < getNIntegrationIntervals(); ++integrationInterval) {
    const unsigned int threadNum = omp_get_thread_num();
    const unsigned int lastRow = mIDCsGrouped.getLastRow();
    unsigned int rowGrouped = 0;
    for (unsigned int iRow = 0; iRow <= lastRow; iRow += mIDCsGrouped.getGroupRows()) {
      // the sectors is divide in to two parts around ylocal=0 to get the same simmetric grouping around ylocal=0
      for (int iYLocalSide = 0; iYLocalSide < 2; ++iYLocalSide) {
        const unsigned int region = mIDCsGrouped.getRegion();
        const unsigned int nPads = Mapper::PADSPERROW[region][iRow] / 2;
        const unsigned int endPads = mIDCsGrouped.getLastPad(iRow) + nPads;

        const unsigned int halfPadsInRow = mIDCsGrouped.getPadsPerRow(rowGrouped) / 2;
        unsigned int padGrouped = iYLocalSide ? halfPadsInRow : halfPadsInRow - 1;
        for (unsigned int ipad = nPads; ipad <= endPads; ipad += mIDCsGrouped.getGroupPads()) {
          const unsigned int endRows = (iRow == lastRow) ? (Mapper::ROWSPERREGION[region] - iRow) : mIDCsGrouped.getGroupRows();
          mRobustAverage[threadNum].clear();
          for (unsigned int iRowMerge = 0; iRowMerge < endRows; ++iRowMerge) {
            const unsigned int iRowTmp = iRow + iRowMerge;
            const auto offs = Mapper::ADDITIONALPADSPERROW[region][iRowTmp] - Mapper::ADDITIONALPADSPERROW[region][iRow];
            const auto padStart = (ipad == 0) ? 0 : offs;
            const unsigned int endPadsTmp = (ipad == endPads) ? (Mapper::PADSPERROW[region][iRowTmp] - ipad) : mIDCsGrouped.getGroupPads() + offs;
            for (unsigned int ipadMerge = padStart; ipadMerge < endPadsTmp; ++ipadMerge) {
              const unsigned int iPadTmp = ipad + ipadMerge;
              const unsigned int iPadSide = iYLocalSide ? iPadTmp : Mapper::PADSPERROW[region][iRowTmp] - iPadTmp - 1;
              const unsigned int indexIDC = integrationInterval * Mapper::PADSPERREGION[region] + Mapper::OFFSETCRULOCAL[region][iRowTmp] + iPadSide;
              mRobustAverage[threadNum].addValue(mIDCsUngrouped[indexIDC] * Mapper::PADAREA[region]);
            }
          }

          switch (paramIDCGroup.Method) {
            case o2::tpc::AveragingMethod::SLOW:
            default:
              mIDCsGrouped(rowGrouped, padGrouped, integrationInterval) = mRobustAverage[threadNum].getFilteredAverage(mSigma);
              break;
            case o2::tpc::AveragingMethod::FAST:
              mIDCsGrouped(rowGrouped, padGrouped, integrationInterval) = mRobustAverage[threadNum].getMean();
              break;
          }

          iYLocalSide ? ++padGrouped : --padGrouped;
        }
      }
      ++rowGrouped;
    }
  }
}

void o2::tpc::IDCAverageGroup::dumpToFile(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "RECREATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}

bool o2::tpc::IDCAverageGroup::setFromFile(const char* fileName, const char* name)
{
  TFile inpf(fileName, "READ");
  IDCAverageGroup* idcAverageGroupTmp{nullptr};
  idcAverageGroupTmp = reinterpret_cast<IDCAverageGroup*>(inpf.GetObjectChecked(name, IDCAverageGroup::Class()));

  if (!idcAverageGroupTmp) {
    LOGP(ERROR, "Failed to load {} from {}", name, inpf.GetName());
    return false;
  }
  setIDCs(idcAverageGroupTmp->getIDCsUngrouped());

  delete idcAverageGroupTmp;
  return true;
}

void o2::tpc::IDCAverageGroup::drawUngroupedIDCs(const unsigned int integrationInterval, const std::string filename) const
{
  const auto coords = o2::tpc::painter::getPadCoordinatesSector();
  TH2Poly* poly = o2::tpc::painter::makeSectorHist("hSector", "Sector;local #it{x} (cm);local #it{y} (cm); #it{IDC}");
  poly->SetContour(255);
  poly->SetTitle(nullptr);
  poly->GetYaxis()->SetTickSize(0.002f);
  poly->GetYaxis()->SetTitleOffset(0.7f);
  poly->GetZaxis()->SetTitleOffset(1.3f);
  poly->SetStats(0);

  TCanvas* can = new TCanvas("can", "can", 2000, 1400);
  can->SetRightMargin(0.14f);
  can->SetLeftMargin(0.06f);
  can->SetTopMargin(0.04f);

  TLatex lat;
  lat.SetTextFont(63);
  lat.SetTextSize(2);

  poly->Draw("colz");
  const unsigned int region = mIDCsGrouped.getRegion();
  for (unsigned int irow = 0; irow < Mapper::ROWSPERREGION[region]; ++irow) {
    for (unsigned int ipad = 0; ipad < Mapper::PADSPERROW[region][irow]; ++ipad) {
      const auto padNum = Mapper::getGlobalPadNumber(irow, ipad, region);
      const auto coordinate = coords[padNum];
      const float yPos = -0.5 * (coordinate.yVals[0] + coordinate.yVals[2]); // local coordinate system is mirrored
      const float xPos = 0.5 * (coordinate.xVals[0] + coordinate.xVals[2]);
      const unsigned int indexIDC = integrationInterval * Mapper::PADSPERREGION[region] + Mapper::OFFSETCRULOCAL[region][irow] + ipad;
      const float idc = mIDCsUngrouped[indexIDC] * Mapper::PADAREA[region];
      poly->Fill(xPos, yPos, idc);
      lat.SetTextAlign(12);
      lat.DrawLatex(xPos, yPos, Form("%i", ipad));
    }
  }

  if (!filename.empty()) {
    can->SaveAs(filename.data());
    delete poly;
    delete can;
  }
}

/// for debugging: creating debug tree for integrated IDCs
/// \param nameFile name of the output file
void o2::tpc::IDCAverageGroup::createDebugTree(const char* nameFile) const
{
  o2::utils::TreeStreamRedirector pcstream(nameFile, "RECREATE");
  pcstream.GetFile()->cd();
  createDebugTree(*this, pcstream);
  pcstream.Close();
}

void o2::tpc::IDCAverageGroup::createDebugTreeForAllCRUs(const char* nameFile, const char* filename)
{
  o2::utils::TreeStreamRedirector pcstream(nameFile, "RECREATE");
  pcstream.GetFile()->cd();
  TFile fInp(filename, "READ");

  for (TObject* keyAsObj : *fInp.GetListOfKeys()) {
    const auto key = dynamic_cast<TKey*>(keyAsObj);
    LOGP(info, "Key name: {} Type: {}", key->GetName(), key->GetClassName());

    if (std::strcmp(o2::tpc::IDCAverageGroup::Class()->GetName(), key->GetClassName()) != 0) {
      LOGP(info, "skipping object. wrong class.");
      continue;
    }

    IDCAverageGroup* idcavg = (IDCAverageGroup*)fInp.Get(key->GetName());
    createDebugTree(*idcavg, pcstream);
    delete idcavg;
  }
  pcstream.Close();
}

void o2::tpc::IDCAverageGroup::createDebugTree(const IDCAverageGroup& idcavg, o2::utils::TreeStreamRedirector& pcstream)
{
  const Mapper& mapper = Mapper::instance();
  unsigned int sector = idcavg.getSector();
  unsigned int cru = sector * Mapper::NREGIONS + idcavg.getRegion();
  const o2::tpc::CRU cruTmp(cru);
  unsigned int region = cruTmp.region();

  for (unsigned int integrationInterval = 0; integrationInterval < idcavg.getNIntegrationIntervals(); ++integrationInterval) {
    const unsigned long padsPerCRU = Mapper::PADSPERREGION[region];
    std::vector<unsigned int> vRow(padsPerCRU);
    std::vector<unsigned int> vPad(padsPerCRU);
    std::vector<float> vXPos(padsPerCRU);
    std::vector<float> vYPos(padsPerCRU);
    std::vector<float> vGlobalXPos(padsPerCRU);
    std::vector<float> vGlobalYPos(padsPerCRU);
    std::vector<float> idcsPerIntegrationInterval(padsPerCRU);        // idcs for one time bin
    std::vector<float> groupedidcsPerIntegrationInterval(padsPerCRU); // idcs for one time bin
    std::vector<float> invPadArea(padsPerCRU);

    for (unsigned int iPad = 0; iPad < padsPerCRU; ++iPad) {
      const GlobalPadNumber globalNum = Mapper::GLOBALPADOFFSET[region] + iPad;
      const auto& padPosLocal = mapper.padPos(globalNum);
      vRow[iPad] = padPosLocal.getRow();
      vPad[iPad] = padPosLocal.getPad();
      vXPos[iPad] = mapper.getPadCentre(padPosLocal).X();
      vYPos[iPad] = mapper.getPadCentre(padPosLocal).Y();
      invPadArea[iPad] = Mapper::PADAREA[region];
      const GlobalPosition2D globalPos = mapper.LocalToGlobal(LocalPosition2D(vXPos[iPad], vYPos[iPad]), cruTmp.sector());
      vGlobalXPos[iPad] = globalPos.X();
      vGlobalYPos[iPad] = globalPos.Y();
      idcsPerIntegrationInterval[iPad] = idcavg.getUngroupedIDCVal(iPad, integrationInterval);
      groupedidcsPerIntegrationInterval[iPad] = idcavg.getGroupedIDCValGlobal(vRow[iPad], vPad[iPad], integrationInterval);
    }

    pcstream << "tree"
             << "cru=" << cru
             << "sector=" << sector
             << "region=" << region
             << "integrationInterval=" << integrationInterval
             << "IDCUngrouped.=" << idcsPerIntegrationInterval
             << "IDCGrouped.=" << groupedidcsPerIntegrationInterval
             << "invPadArea.=" << invPadArea
             << "pad.=" << vPad
             << "row.=" << vRow
             << "lx.=" << vXPos
             << "ly.=" << vYPos
             << "gx.=" << vGlobalXPos
             << "gy.=" << vGlobalYPos
             << "\n";
  }
}

void o2::tpc::IDCAverageGroup::setIDCs(const std::vector<float>& idcs)
{
  mIDCsUngrouped = idcs;
  mIDCsGrouped.resize(getNIntegrationIntervals());
}

void o2::tpc::IDCAverageGroup::setIDCs(std::vector<float>&& idcs)
{
  mIDCsUngrouped = std::move(idcs);
  mIDCsGrouped.resize(getNIntegrationIntervals());
}

unsigned int o2::tpc::IDCAverageGroup::getNIntegrationIntervals() const
{
  return mIDCsUngrouped.size() / Mapper::PADSPERREGION[mIDCsGrouped.getRegion()];
}

float o2::tpc::IDCAverageGroup::getUngroupedIDCVal(const unsigned int localPadNumber, const unsigned int integrationInterval) const
{
  return mIDCsUngrouped[localPadNumber + integrationInterval * Mapper::PADSPERREGION[mIDCsGrouped.getRegion()]];
}

unsigned int o2::tpc::IDCAverageGroup::getUngroupedIndex(const unsigned int ulrow, const unsigned int upad, const unsigned int integrationInterval) const
{
  return integrationInterval * Mapper::PADSPERREGION[mIDCsGrouped.getRegion()] + Mapper::OFFSETCRULOCAL[mIDCsGrouped.getRegion()][ulrow] + upad;
}

unsigned int o2::tpc::IDCAverageGroup::getUngroupedIndexGlobal(const unsigned int ugrow, const unsigned int upad, const unsigned int integrationInterval) const
{
  return integrationInterval * Mapper::PADSPERREGION[mIDCsGrouped.getRegion()] + Mapper::OFFSETCRUGLOBAL[ugrow] + upad;
}
