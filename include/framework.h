#ifndef FRAMEWORK_CLASS_H
#define FRAMEWORK_CLASS_H

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <random>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <getopt.h>
#include <chrono>
#include <iomanip>
#include <ctime>


class framework{
public:

  struct struct_Sys_properties{ /*! */
    float x,                    /*! */
          y,                    /*! */
          z,                    /*! */
          R;                    /*! */

    bool tSpin;                 /*! */
    int spin;
    std::vector<bool> sBorder;
  };

  struct Fs_general_info{    /*! structure to contain global properties or magnitude of the system. */
    double FsTemperature,    /*! system temperature */
           FsPressure,       /*! system pressure */
           FsVolume,         /*! system volume */
           FsL,              /*! L represent the length of a system side */
           FsTotalEnergy,      /*! Absolute value of total energy of the system */
           FsTotalMMoment,     /*! Absolute value of total magnetic moment system, represent the summation of al n = 1 or 0 spin states  */
           Fsk_nn,             /*! Fsk_nn represent for the moment the spring constant, but will be a function as soon as possible, is in the order of K_B (e-23) */
           FsMinT,             /*! min value of temperature cycle*/
           FsMaxT,             /*! max value of temperature cycle */
           FsTsteps,
           FsD,                /*! D represent the value of energy between HS-LS state, is in the order of K_B (e-23) */
           Fsj_isin,           /*! Fsj_isin for this moment represent the value interaction in on a system isin-like mode, is in the order of K_B (e-23) */
           FsTEnergy_error,    /*! FsTEnergy_error is used as the representation absolute energy error per temperature   */
           FsTMMoment_error,   /*! FsTMMoment_error is used as the representation absolute magnetic moment error per temperature   */
           FsR_hs,             /*! usually 0.418*/
           FsR_ls,
           Fsdelta;            /*! initial distance between atoms(default=0.9). */
           std::vector<double> FsLlist;

    int FsTotalAtoms,       /*! total of atoms set, will be FsNi^3, this variable is used like a N (total atoms) */
        FsNi,               /*! number of atoms per side sqrt(N)*/
        Fscycles,           /*! amount of cycles in the Monte Carlos method  */
        FsScreenrate,
        FsFilerate,
        FsSlice;

  };

  struct statistic_s{
    std::vector<double> sSEnergy_vector;         /*! FsEnergy_vector, vector that will contain each -Fscycles- energy values to do statistic operations */
    std::vector<double> sSLiveAverageEnergy;
    std::vector<double> sSEnergy_STD_error_vector;
    double sSLiveAveE;
    double sSFinalAverageEnergy;
    double sSEnergy_STD_error;

    std::vector<double> sSMMoment_vector;        /*! FsMMoment_vector, vector that will contain each -Fscycles- magnetic moment values to do statistic operations */
    std::vector<double> sSLiveAverageMMoment;
    std::vector<double> sSMoment_STD_error_vector;
    double sSLiveAveMM;
    double sSFinalAverageMMoment;
    double sSMoment_STD_error;
  };

  struct s_options{
    int a_temperatureCycle,
        b_temperatureCycle,
        totalEnergy_totalMoment,
        routine,
        metropolisMethod,
        swithEnergy,
        print_fileout,
        det_VandL,
        print_out,
        liveAverage,
        xyz_print,
        temperature_file_print,
        cycle_file_print,
        test_function,
        about_system_inf,
        VMD_EXPORT_F,
        open_parammeters_file;
  };

  float R_hs = 0.418;     /*! R_ls=0.38 definition, with the relation R_hs/R_ls=1.1. Take the distance between couples of atoms as 1, x = 0.202= 1 - R_hs + R_hs*/
  float L =0;             /*! variable that will be remove*/

  framework();
  framework(std::vector<struct_Sys_properties> &,    /*!vector type structure struct_Sys_properties*/
            std::vector<std::vector<int> > &,        /*!2D int vector */
            const double,                            /*const double Temperature*/
            const double,                            /*const double Pressure*/
            const int);                              /*!Atoms, where the total are N*N*N.*/
  ~framework();

  int f_define_system_info(Fs_general_info &,
                           statistic_s &,
                           s_options &);

  int f_write_vector (std::vector<struct_Sys_properties> & ,        /*!vector type structure*/
                      const int);                                   /*!Atoms, where the total are N*N*N.*/

  int f_read_vector (const std::vector<struct_Sys_properties> &,    /*!vector type structure*/
                     const int);                                    /*!int.Atoms, where the total are N*N*N.*/

  int f_print_fileout (const std::vector<struct_Sys_properties> &,  /*!vector type structure*/
                       const std::vector<std::vector<int> > &,      /*!2D int vector */
                       Fs_general_info &,                           /*!*/
                       const statistic_s & ,
                       const int,                                   /*!int.option to switch */
                       const int);                                  /*!*/

  int f_metropolisMethod (std::vector<struct_Sys_properties> &,              /*!vector type structure*/
                          const std::vector<std::vector<int> > &,            /*!2D int vector */
                          Fs_general_info &,
                          statistic_s &,
                          s_options &,
                          const int);                               /*!*/

  double f_switchEnergy(const std::vector<struct_Sys_properties> &,           /*!const.vector type structure*/
                       const std::vector<std::vector<int> > &,               /*!2D int vector */
                       struct_Sys_properties&,                               /*!framework::struct_Sys_properties temporal recipient*/
                       Fs_general_info &,                                    /*!*/
                       const int,
                       const int);                                           /*!cosnt.int.Index to do de switch*/

  double f_cell_volume(const std::vector<struct_Sys_properties> &);      /*!const. vector type structure*/

  double f_det_VandL (const std::vector<struct_Sys_properties> &,        /*!const.vector type structure*/
                      Fs_general_info &,                                 /*!*/
                      const int);                                        /*!*/

  int f_routine (std::vector<struct_Sys_properties> &,             /*!vector type structure*/
                 const std::vector<std::vector<int> > &,           /*!2D int vector */
                 Fs_general_info &,                                /*! system info*/
                 statistic_s & ,
                 s_options &,
                 const int,                                        /*!const int option */
                 const int);                                        /*!const int number of cycles*/

  int create_lattice (std::vector<struct_Sys_properties> &,        /*!vector type structure*/
                 std::vector<std::vector<int> > &,                 /*!2D int vector */
                 framework::Fs_general_info &,                     /*! */
                 const float);                                     /*!const.float.R_hs, high spin Radio*/


  int f_define_system_info (Fs_general_info & system_info,
                            statistic_s &,
                            s_options & ,
                            const double,
                            const double,
                            const double,
                            const double,
                            const double,
                            const double,
                            const double,
                            const double,
                            const double,
                            const int,
                            const int);                                 /*!*/

  int f_define_system_info (Fs_general_info &,
                            statistic_s &,
                            s_options &,
                            const double,
                            const int,
                            const int);

  int f_temperatureCycle (std::vector<struct_Sys_properties> &,                /*!vector type structure*/
                          std::vector<std::vector<int> > &,                    /*!2D int vector */
                          Fs_general_info&,                                    /*!*/
                          statistic_s &,
                          s_options &,
                          const double,                                        /*!*/
                          const double,                                        /*!*/
                          const double);                                       /*!*/
  int f_temperatureCycle(std::vector<struct_Sys_properties> &,
                         std::vector<std::vector<int> > &,
                         Fs_general_info&,
                         statistic_s &,
                         s_options &,
                         const int);

  double f_totalEnergy_totalMoment (std::vector<struct_Sys_properties> &,     /*!vector type structure*/
                                 const std::vector<std::vector<int> > &,                /*!2D int vector */
                                 Fs_general_info&,                     /*!*/
                                 statistic_s &,
                                 const int,
                                 const int /*!opt*/ );                            /*!*/

  int f_print_out(const std::vector<struct_Sys_properties> &,              /*!vector type structure*/
                  const std::vector<std::vector<int> > &,                  /*!2D int vector */
                  const Fs_general_info &,                                 /*!opt*/
                  const int);                                              /*!opt*/

  int f_live_statistic (Fs_general_info&,
                        statistic_s &);

  double f_liveAverage (statistic_s &,
                        const int,
                        const int);

  int f_detect_border(std::vector<struct_Sys_properties> &,
                      const framework::Fs_general_info &);

  int f_xyz_print(const std::vector<struct_Sys_properties> &,
                  const Fs_general_info &,
                  const statistic_s &,
                  const int,
                  const int);

  int f_temperature_file_print(const std::vector<struct_Sys_properties> &,
                  const Fs_general_info &,
                  const statistic_s &,
                  const int);

  int f_cycle_file_print(const std::vector<struct_Sys_properties> &,
                         const Fs_general_info &,
                         const statistic_s &,
                         const s_options &,
                         const int,
                         const int);

  double f_tow_p_distance(const std::vector<struct_Sys_properties> &,
                          const Fs_general_info &,
                          const std::vector<std::vector<int> > &,
                          const int,
                          const int);

  double f_tow_p_distance(const std::vector<struct_Sys_properties> &,
                          const struct_Sys_properties &,
                          const Fs_general_info &,
                          const std::vector<std::vector<int> > &,
                          const int,
                          const int);

  double f_test_function(std::vector<struct_Sys_properties> &,           /*!const.vector type structure*/
                         Fs_general_info &,
                         s_options & system_options,                                  /*!*/
                         const std::vector<std::vector<int> > &,               /*!2D int vector */
                         statistic_s &,
                         const int,                              /*!framework::struct_Sys_properties temporal recipient*/
                         const int);

  void f_about_system_inf(const Fs_general_info &,
                          const int);

  int f_option_selectior(const Fs_general_info &);

  int f_creating_headfileout(const int opt);

  void f_remove_olf_files();

  double f_STD_desviation (Fs_general_info&,
                           statistic_s &,
                           const int,
                           const int opt);










private:
  std::vector<struct_Sys_properties> main_system;                          /*! */
  std::vector<std::vector<int> > main_system_neighbour;                    /*! */
  Fs_general_info main_system_info;                                             /*! */
  statistic_s main_statistic_info;
  s_options main_system_options;
};


/******************************************/
/*  PROTOTYPE AREA OF AUXILIAR FUNCTIONS  */
/*                                        */
/******************************************/

int int_randomNumber_generator(const int,                                     /*!int.min*/
                               const int);                                    /*!int.max*/

double real_randomNumber_generator(const int,                                 /*!int.min*/
                                   const int);                                /*!int.max*/

int randomSamplePosition(std::vector<int> &,                                  /*!vector int*/
                         const int);                                          /*!const int.N*/

int neighbourDefine (const std::vector<framework::struct_Sys_properties> &,   /*!const. vector type structure*/
                     std::vector<std::vector<int> > &,                        /*!2D int vector */
                     const int);                                              /*!const.int.N*/

int corner_f_lattice (const std::vector<framework::struct_Sys_properties> &,  /*!const. vector type structure*/
                      std::vector<int> &,                                     /*int. vector to put the corner values of lattice*/
                      const int);                                             /*!const.int.N*/

int VMD_EXPORT_F (const std::vector<framework::struct_Sys_properties> &,      /*! */
                  const int,                                                  /*! */
                  const int);                                                 /*! Choice options: 1 to remove old files or 2 to create the vmd files. */
int parse_commandline(int argc,
                      char **argv,
                      framework::Fs_general_info &,
                      framework::s_options & );
void help ();

#endif // MSYSTEM_CLASS_H
