#ifndef DATA_STRUCT_H_
#define DATA_STRUCT_H_

#include <cstdlib>
#include <iostream>
#include <string>
#include <array>
#include <vector>

namespace iSS_data {
    const double hbarC=0.197327053;  //GeV*fm

    typedef std::array<float, 4> Vec4;
    typedef std::array<double, 4> VecD4;
    typedef std::array<double, 10> ViscousVec;

    const int NUMBER_OF_LINES_TO_WRITE = 100000;
    const int AMOUNT_OF_OUTPUT = 0;
}


enum class AfterburnerType {
    UrQMD,
    SMASH,
    JAM,
    PDG_Decay,
};


//! The particle list's name in the table files iSS reads:
//! iSS_tables/pdg-<name>.dat (species and decays),
//! iSS_tables/chosen_particles_<name>.dat (species sampled) and
//! iSS_tables/EOS_tables/HRGEOS_PST-<name>.dat (its hadron gas, regulateEOS).
//! By default (parameter particle_list = 0) the list of the afterburner;
//! particle_list = 21 selects PDG2021+ (arXiv:2309.01737) whatever the
//! afterburner, whose own list must then be the same.  chosen_list asks for
//! the name of the chosen-particle file, which differs for s95p-v1.
inline std::string particle_list_name(int particle_list,
                                      AfterburnerType afterburner,
                                      bool chosen_list = false) {
    if (particle_list == 21) return "PDG2021Plus";
    if (particle_list != 0) {
        std::cerr << "[Error] unknown particle_list = " << particle_list
                  << " (0: the afterburner's list, 21: PDG2021Plus)"
                  << std::endl;
        exit(1);
    }
    if (afterburner == AfterburnerType::SMASH) return "SMASH";
    if (afterburner == AfterburnerType::UrQMD) return "urqmd_v3.3+";
    return chosen_list ? "s95p-v1" : "s95pv1";
}


typedef struct {
    int decay_Npart;
    double branching_ratio;
    int decay_part[5];
} decay_channel_info;


typedef struct {
    int monval;     // Monte Carlo number according PDG
    std::string name;
    double mass;
    double width;
    int gspin;      // spin degeneracy
    int baryon;
    int strange;
    int charm;
    int bottom;
    int gisospin;   // isospin degeneracy
    int charge;
    int decays;     // amount of decays listed for this resonance
    int stable;     // defines whether this particle is considered as stable
    std::vector<decay_channel_info*> decay_channels;
    int sign;                   // Bose-Einstein or Dirac-Fermi statistics
} particle_info;


typedef struct {
    float tau, xpt, ypt, eta;
    float cosh_eta, sinh_eta;
    float da0, da1, da2, da3;
    float u0, u1, u2, u3;
    float Edec, Tdec, Pdec;
    float Bn, Qn, Sn, muB, muS, muQ;
    float pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33;
    float bulkPi;
    float qmu0, qmu1, qmu2, qmu3;
    std::vector<float> particle_mu_PCE;
    std::vector<float> vorticity_arr;
    std::vector<double> visCoeffs;
} FO_surf;


typedef struct {
    float tau, xpt, ypt, eta;
    iSS_data::Vec4 da_mu_LRF;
    iSS_data::Vec4 u_tz;
    float Edec, Tdec, Pdec;
    float Bn, Qn, Sn, muB, muS, muQ;
    float bulkPi;
    float piLRF_xx, piLRF_xy, piLRF_xz, piLRF_yy, piLRF_yz;
    float qmuLRF_x, qmuLRF_y, qmuLRF_z;
    std::vector<float> particle_mu_PCE;
    std::vector<double> visCoeffs;
} FO_surf_LRF;


struct iSS_Hadron {
     int pid;
     float mass;
     float E, px, py, pz;
     float t, x, y, z;
};

#endif  // DATA_STRUCT_H_
