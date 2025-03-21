#ifndef FT3MODULE_H
#define FT3MODULE_H

#include <TGeoVolume.h>
#include <string>

class FT3Module {
    
public:

    static void initialize_materials();
    static TGeoMaterial* siliconMat;
    static TGeoMedium* siliconMed;
    static TGeoMaterial* copperMat;
    static TGeoMedium* copperMed;
    static TGeoMaterial* kaptonMat;
    static TGeoMedium* kaptonMed;
    static TGeoMaterial* epoxyMat;
    static TGeoMedium* epoxyMed;
    static TGeoMaterial* AluminumMat;
    static TGeoMedium* AluminumMed; 


    const char* mDetName; 

    static void createModule(double mZ, int layerNumber, int direction, double Rin, double Rout, double overlap, const std::string& face, const std::string& layout_type, TGeoVolume* motherVolume);


private:

    static void create_layout(double mZ, int layerNumber, int direction, double Rin, double Rout, double overlap, const std::string& face, const std::string& layout_type, TGeoVolume* motherVolume);

};

#endif // FT3MODULE_H
