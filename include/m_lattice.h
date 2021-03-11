#ifndef m_lattice_H
#define m_lattice_H


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




class m_lattice
{
  public:

    m_lattice();
    virtual ~m_lattice();

    struct s_three_D_system{ /*! */
      float x,                    /*! */
            y,                    /*! */
            z,                    /*! */
            R;                    /*! */

      bool tSpin;                 /*! */
      int spin;
      std::vector<bool> sBorder;
    };

    struct s_two_D_system{ /*! */
      float x,                    /*! */
            y,                    /*! */
            R;                    /*! */
      bool tSpin;                 /*! */
      int spin;
      std::vector<bool> sBorder;
      std::vector<int> two_D_neighbour;
    };

    struct s_one_D_system{ /*! */
      float x,                    /*! */
            R;                    /*! */
      bool tSpin;                 /*! */
      int spin;
      std::vector<bool> sBorder;
    };

    struct s_general_info{    /*! structure to contain global properties or magnitude of the system. */
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

    int run(std::vector<s_one_D_system> &, s_general_info &);

    int f_create_3D_lattice();
    int f_create_2D_lattice(std::vector<s_two_D_system> &, s_general_info &);
    int f_create_1D_lattice(std::vector<s_one_D_system> &, s_general_info &);

  protected:

  private:
};

#endif // m_lattice



