#include<framework.h>

/******************/
/*   CLASS AREA   */
/******************/


/************************************/
/* CONSTRUCTOR AND DESTRUCTOR AREA  */
/************************************/

framework::framework()
{
  /*!---------------- initialising system info ----------------*/
  main_system_info.FsTemperature = 0.1;  /*! system temperature. */
  main_system_info.FsPressure = 0.01;    /*! system pressure. */
  main_system_info.FsVolume = 0;         /*! system volume. */
  main_system_info.FsL = 0;              /*! FsL will be modified or calculate. */
  main_system_info.FsTotalEnergy = 0;    /*! will be modified or calculate. */
  main_system_info.FsTotalMMoment = 0;   /*! will be modified or calculate. */
  main_system_info.Fsk_nn = 1;           /*! Fsk_nn begin in 1 but will be change, but will be a function as soon as possible, is in the order of K_B (e-23). */
  main_system_info.FsMinT = 1.0;         /*! min value of temperature cycle. */
  main_system_info.FsMaxT = 10.0;         /*! max value of temperature cycle. */
  main_system_info.FsTsteps = 0.2;
  main_system_info.FsD = 1;              /*! D represent the value of energy between HS-LS state, is in the order of K_B (e-23). */
  main_system_info.Fsj_isin = 1;           /*! Fsj_isin for this moment represent the value interaction in on a system isin-like mode, is in the order of K_B (e-23). */
  main_system_info.FsTEnergy_error = 0;    /*! FsTEnergy_error is used as the representation absolute energy error per temperature. */
  main_system_info.FsTMMoment_error = 0;   /*! FsTMMoment_error is used as the representation absolute magnetic moment error per temperature. */
  main_system_info.FsR_hs = 0.418;         /*! with the relation R_hs/R_ls=1.1. Take the distance between couples of atoms as 1, x = 0.202= 1 - R_hs + R_hs.*/
  main_system_info.FsR_ls = 0.38;          /*! radio to correspond to low-spin state. */
  main_system_info.FsSlice = 500;

  main_system_info.FsTotalAtoms = 27;     /*! total of atoms set, will be FsNi^3, this variable is used like a N (total atoms) */
  main_system_info.FsNi = 3;               /*! number of atoms per side sqrt(N)*/
  main_system_info.Fscycles = 1000;         /*! amount of cycles in the Monte Carlos method  */
  main_system_info.Fsdelta = 0.9;

  main_statistic_info.sSEnergy_vector.resize(main_system_info.Fscycles);    /*! sSEnergy_vector, vector that will contain each -Fscycles- energy values to do statistic operations */
  main_statistic_info.sSLiveAveE = 0;
  main_statistic_info.sSFinalAverageEnergy = 0;
  main_statistic_info.sSMMoment_vector.resize(main_system_info.Fscycles);  /*! sSMMoment_vector, vector that will contain each -Fscycles- magnetic moment values to do statistic operations */
  main_statistic_info.sSLiveAveMM = 0;
  main_statistic_info.sSFinalAverageMMoment = 0;

  main_statistic_info.sSEnergy_vector.resize(main_system_info.Fscycles);                  //!ver si esta inicialization se hace antes
  std::fill(main_statistic_info.sSEnergy_vector.begin(), main_statistic_info.sSEnergy_vector.end(), 0);
  main_statistic_info.sSLiveAverageEnergy.resize(main_system_info.Fscycles);              //!
  std::fill(main_statistic_info.sSLiveAverageEnergy.begin(), main_statistic_info.sSLiveAverageEnergy.end(), 0);

  main_statistic_info.sSMMoment_vector.resize(main_system_info.Fscycles);                 //!ver si esta inicialization se hace antes
  std::fill(main_statistic_info.sSMMoment_vector.begin(), main_statistic_info.sSMMoment_vector.end(), 0);
  main_statistic_info.sSLiveAverageMMoment.resize(main_system_info.Fscycles);
  std::fill(main_statistic_info.sSLiveAverageMMoment.begin(), main_statistic_info.sSLiveAverageMMoment.end(), 0);

  main_statistic_info.sSEnergy_STD_error_vector.resize(main_system_info.Fscycles);
  std::fill(main_statistic_info.sSEnergy_STD_error_vector.begin(), main_statistic_info.sSEnergy_STD_error_vector.end(), 0);
  main_statistic_info.sSMoment_STD_error_vector.resize(main_system_info.Fscycles);
  std::fill(main_statistic_info.sSMoment_STD_error_vector.begin(), main_statistic_info.sSMoment_STD_error_vector.end(), 0);

  main_system_options.open_parammeters_file = 1;

 /*!---------------- initialising system system and neighbour ----------------*/
  main_system.resize(main_system_info.FsTotalAtoms);

  main_system_neighbour.resize(main_system_info.FsTotalAtoms);

  int flag = 0;
  for (int z = 0; z < main_system_info.FsNi; ++z){
    for (int y = 0; y < main_system_info.FsNi; ++y){
      for (int x = 0; x < main_system_info.FsNi; ++x){
       main_system[flag].x = (2 * x + 1) * (main_system_info.Fsdelta/2 + main_system_info.FsR_hs); //
       main_system[flag].y = (2 * y + 1) * (main_system_info.Fsdelta/2 + main_system_info.FsR_hs);
       main_system[flag].z = (2 * z + 1) * (main_system_info.Fsdelta/2 + main_system_info.FsR_hs);
       main_system[flag].tSpin=1;
       main_system[flag].spin=1;
       main_system[flag].R = main_system_info.FsR_hs;
       main_system[flag].sBorder.resize(6);
       main_system_neighbour[flag].resize(6);

       //std::cout<<flag<<std::endl;
       flag++;
      }
    }
  }
  // stimating L and the volumen.
  //main_system_info.FsL = f_det_VandL (main_system, main_system_info, 2);
  main_system_info.FsL = 0.9 * main_system_info.FsNi;
  main_system_info.FsVolume = main_system_info.FsL * main_system_info.FsL * main_system_info.FsL;
}

framework::framework(std::vector<struct_Sys_properties> & system_in,
	                   std::vector<std::vector<int> > &neighbour_vector,
	                   const double Temperature,
	                   const double Pressure,
	                   const int n)
{
 /*!---------------- initialising system info ----------------*/
  main_system_info.FsTemperature = Temperature;  /*! system temperature. */
  main_system_info.FsPressure = Pressure;        /*! system pressure. */
  main_system_info.FsVolume = 0;                 /*! system volume. */
  main_system_info.FsL = 0;                      /*! FsL will be modified or calculate. */
  main_system_info.FsTotalEnergy = 0;            /*! will be modified or calculate. */
  main_system_info.FsTotalMMoment = 0;           /*! will be modified or calculate. */
  main_system_info.Fsk_nn = 1;                   /*! Fsk_nn begin in 1 but will be change, but will be a function as soon as possible, is in the order of K_B (e-23). */
  main_system_info.FsMinT = 1.0;                 /*! min value of temperature cycle. */
  main_system_info.FsMaxT = 10.0;                 /*! max value of temperature cycle. */
  main_system_info.FsTsteps = 0.2;
  main_system_info.FsD = 1;                      /*! D represent the value of energy between HS-LS state, is in the order of K_B (e-23). */
  main_system_info.Fsj_isin = 1;           /*! Fsj_isin for this moment represent the value interaction in on a system isin-like mode, is in the order of K_B (e-23). */
  main_system_info.FsTEnergy_error = 0;    /*! FsTEnergy_error is used as the representation absolute energy error per temperature. */
  main_system_info.FsTMMoment_error = 0;   /*! FsTMMoment_error is used as the representation absolute magnetic moment error per temperature. */
  main_system_info.FsR_hs = 0.418;         /*! with the relation R_hs/R_ls=1.1. Take the distance between couples of atoms as 1, x = 0.202= 1 - R_hs + R_hs.*/
  main_system_info.FsR_ls = 0.38;          /*! radio to correspond to low-spin state. */

  main_system_info.FsTotalAtoms = n * n * n;     /*! total of atoms set, will be FsNi^3, this variable is used like a N (total atoms) */
  main_system_info.FsNi = n;                     /*! number of atoms per side sqrt(N)*/
  main_system_info.Fscycles = 1000;              /*! amount of cycles in the Monte Carlos method  */
  main_system_info.FsSlice = 500;
  main_system_info.Fsdelta = 0.9;

  main_statistic_info.sSEnergy_vector.resize(main_system_info.Fscycles);    /*! sSEnergy_vector, vector that will contain each -Fscycles- energy values to do statistic operations */
  main_statistic_info.sSLiveAveE = 0;
  main_statistic_info.sSFinalAverageEnergy = 0;
  main_statistic_info.sSMMoment_vector.resize(main_system_info.Fscycles);  /*! sSMMoment_vector, vector that will contain each -Fscycles- magnetic moment values to do statistic operations */
  main_statistic_info.sSLiveAveMM = 0;
  main_statistic_info.sSFinalAverageMMoment = 0;

  main_statistic_info.sSEnergy_vector.resize(main_system_info.Fscycles);                  //!ver si esta inicialization se hace antes
  std::fill(main_statistic_info.sSEnergy_vector.begin(), main_statistic_info.sSEnergy_vector.end(), 0);
  main_statistic_info.sSLiveAverageEnergy.resize(main_system_info.Fscycles);              //!
  std::fill(main_statistic_info.sSLiveAverageEnergy.begin(), main_statistic_info.sSLiveAverageEnergy.end(), 0);

  main_statistic_info.sSMMoment_vector.resize(main_system_info.Fscycles);                 //!ver si esta inicialization se hace antes
  std::fill(main_statistic_info.sSMMoment_vector.begin(), main_statistic_info.sSMMoment_vector.end(), 0);
  main_statistic_info.sSLiveAverageMMoment.resize(main_system_info.Fscycles);
  std::fill(main_statistic_info.sSLiveAverageMMoment.begin(), main_statistic_info.sSLiveAverageMMoment.end(), 0);

  main_statistic_info.sSEnergy_STD_error_vector.resize(main_system_info.Fscycles);
  std::fill(main_statistic_info.sSEnergy_STD_error_vector.begin(), main_statistic_info.sSEnergy_STD_error_vector.end(), 0);
  main_statistic_info.sSMoment_STD_error_vector.resize(main_system_info.Fscycles);
  std::fill(main_statistic_info.sSMoment_STD_error_vector.begin(), main_statistic_info.sSMoment_STD_error_vector.end(), 0);


  main_system_options.open_parammeters_file = 0;

  /*!---------------- creating system system and neighbour (class and external) ----------------*/
  main_system.resize(main_system_info.FsTotalAtoms);
  main_system_neighbour.resize(main_system_info.FsTotalAtoms);

  system_in.resize(main_system_info.FsTotalAtoms);           /*! creating the external vector system.*/
  neighbour_vector.resize(main_system_info.FsTotalAtoms);    /*! creating the external neighbour vector system.*/


  /*!---------------- initialising system system and neighbour ----------------*/
  int flag = 0;
  for (int z = 0; z < main_system_info.FsNi; ++z){
    for (int y = 0; y < main_system_info.FsNi; ++y){
      for (int x = 0; x < main_system_info.FsNi; ++x){
       main_system[flag].x = (2 * x + 1) * (0.202/2 + main_system_info.FsR_hs);
       main_system[flag].y = (2 * y + 1) * (0.202/2 + main_system_info.FsR_hs);
       main_system[flag].z = (2 * z + 1) * (0.202/2 + main_system_info.FsR_hs);
       main_system[flag].tSpin=1;
       main_system[flag].spin=1;
       main_system[flag].sBorder.resize(6);
       if (main_system[flag].tSpin == 1) main_system[flag].R = main_system_info.FsR_hs;
         else  main_system[flag].R = main_system_info.FsR_hs/1.1;
       main_system_neighbour[flag].resize(6);
       //std::cout<<flag<<std::endl;
       flag++;
      }
    }
  }

  // stimating L and the volumen.
  //main_system_info.FsL = f_det_VandL (main_system, main_system_info, 2);
  main_system_info.FsL = 0.9 * main_system_info.FsNi;
  main_system_info.FsVolume = main_system_info.FsL * main_system_info.FsL * main_system_info.FsL;

}

framework::~framework(){}


/*********************/
/*  CLASS FUNCTIONS  */
/*********************/

int framework::f_define_system_info(framework::Fs_general_info & system_info,
                                    statistic_s & statistic_info,
                                    s_options & system_options,
	                                  const double Temperature,
	                                  const double Pressure,
	                                  const double K_nn,
	                                  const double Tmin,
	                                  const double Tmax,
	                                  const double Tsteps ,
	                                  const double D,
	                                  const double j_isin,
	                                  const double R_hs,
	                                  const int cycles,
	                                  const int n)
{
  /*!---------------- initialising system info ----------------*/
  system_info.FsTemperature = Temperature;  /*! system temperature. */
  system_info.FsPressure = Pressure;        /*! system pressure. */
  system_info.FsVolume = 0;                 /*! system volume. */
  system_info.FsL = 0;                      /*! FsL will be modified or calculate. */
  system_info.FsTotalEnergy = 0;            /*! will be modified or calculate. */
  system_info.FsTotalMMoment = 0;           /*! will be modified or calculate. */
  system_info.Fsk_nn = K_nn;                /*! Fsk_nn begin in 1 but will be change, but will be a function as soon as possible, is in the order of K_B (e-23). */
  system_info.FsMinT = Tmin;                /*! min value of temperature cycle. */
  system_info.FsMaxT = Tmax;                /*! max value of temperature cycle. */
  system_info.FsTsteps = Tsteps;
  system_info.FsD = D;                      /*! D represent the value of energy between HS-LS state, is in the order of K_B (e-23). */
  system_info.Fsj_isin = j_isin;            /*! Fsj_isin for this moment represent the value interaction in on a system isin-like mode, is in the order of K_B (e-23). */
  system_info.FsTEnergy_error = 0;          /*! FsTEnergy_error is used as the representation absolute energy error per temperature. */
  system_info.FsTMMoment_error = 0;         /*! FsTMMoment_error is used as the representation absolute magnetic moment error per temperature. */
  system_info.FsR_hs = R_hs;                /*! with the relation R_hs/R_ls=1.1. Take the distance between couples of atoms as 1, x = 0.202= 1 - R_hs + R_hs.*/
  system_info.FsR_ls = R_hs/1.1;            /*! radio to correspond to low-spin state. */
  system_info.Fsdelta = 0.9;

  system_info.FsTotalAtoms = n * n * n;     /*! total of atoms set, will be FsNi^3, this variable is used like a N (total atoms) */
  system_info.FsNi = n;                     /*! number of atoms per side sqrt(N)*/
  system_info.Fscycles = cycles;            /*! amount of cycles in the Monte Carlos method  */
  system_info.FsScreenrate = 500;
  system_info.FsFilerate = 100;
  system_info.FsSlice = 500;

  statistic_info.sSEnergy_vector.resize(system_info.Fscycles);    /*! sSEnergy_vector, vector that will contain each -Fscycles- energy values to do statistic operations */
  statistic_info.sSLiveAveE = 0;
  statistic_info.sSFinalAverageEnergy = 0;
  statistic_info.sSMMoment_vector.resize(system_info.Fscycles);  /*! sSMMoment_vector, vector that will contain each -Fscycles- magnetic moment values to do statistic operations */
  statistic_info.sSLiveAveMM = 0;
  statistic_info.sSFinalAverageMMoment = 0;

  // defining options by default
  system_options.b_temperatureCycle = 1;
  system_options.totalEnergy_totalMoment = 0;
  system_options.routine = 3;
  system_options.metropolisMethod = 4;
  system_options.swithEnergy = 3;
  system_options.print_fileout = 1;
  system_options.det_VandL = 1;
  system_options.print_out = 1;
  system_options.liveAverage = 2;
  system_options.xyz_print = 1;
  system_options.temperature_file_print = 1;
  system_options.cycle_file_print = 1;
  system_options.test_function = 2;
  system_options.about_system_inf = 1;
  system_options.VMD_EXPORT_F = 1;

  statistic_info.sSEnergy_vector.resize(system_info.Fscycles);                  //!ver si esta inicialization se hace antes
  std::fill(statistic_info.sSEnergy_vector.begin(), statistic_info.sSEnergy_vector.end(), 0);
  statistic_info.sSLiveAverageEnergy.resize(system_info.Fscycles);              //!
  std::fill(statistic_info.sSLiveAverageEnergy.begin(), statistic_info.sSLiveAverageEnergy.end(), 0);

  statistic_info.sSMMoment_vector.resize(system_info.Fscycles);                 //!ver si esta inicialization se hace antes
  std::fill(statistic_info.sSMMoment_vector.begin(), statistic_info.sSMMoment_vector.end(), 0);
  statistic_info.sSLiveAverageMMoment.resize(system_info.Fscycles);
  std::fill(statistic_info.sSLiveAverageMMoment.begin(), statistic_info.sSLiveAverageMMoment.end(), 0);

  statistic_info.sSEnergy_STD_error_vector.resize(system_info.Fscycles);
  std::fill(statistic_info.sSEnergy_STD_error_vector.begin(), statistic_info.sSEnergy_STD_error_vector.end(), 0);
  statistic_info.sSMoment_STD_error_vector.resize(system_info.Fscycles);
  std::fill(statistic_info.sSMoment_STD_error_vector.begin(), statistic_info.sSMoment_STD_error_vector.end(), 0);

  return 0;
}

int framework::f_define_system_info(framework::Fs_general_info & system_info,
                                    statistic_s & system_statistic,
                                    s_options & system_options,
	                                  const double Temperature,
	                                  const int cycles,
	                                  const int n)
{
  /*!---------------- initialising system info ----------------*/
  system_info.FsTemperature = Temperature;  /*! system temperature. */
  system_info.FsPressure = 0.01;            /*! system pressure. */
  system_info.FsVolume = 0;                 /*! system volume. */
  system_info.FsL = 0;                      /*! FsL will be modified or calculate. */
  system_info.FsTotalEnergy = 0;            /*! will be modified or calculate. */
  system_info.FsTotalMMoment = 0;           /*! will be modified or calculate. */
  system_info.Fsk_nn = 10;                  /*! Fsk_nn begin in 1 but will be change, but will be a function as soon as possible, is in the order of K_B (e-23). */
  system_info.FsMinT = 1.0;                 /*! min value of temperature cycle. */
  system_info.FsMaxT = 10.0;                  /*! max value of temperature cycle. */
  system_info.FsTsteps = 0.2;
  system_info.FsD = 1.0;                    /*! D represent the value of energy between HS-LS state, is in the order of K_B (e-23). */
  system_info.Fsj_isin = 1.0;              /*! Fsj_isin for this moment represent the value interaction in on a system isin-like mode, is in the order of K_B (e-23). */
  system_info.FsTEnergy_error = 0;          /*! FsTEnergy_error is used as the representation absolute energy error per temperature. */
  system_info.FsTMMoment_error = 0;         /*! FsTMMoment_error is used as the representation absolute magnetic moment error per temperature. */
  system_info.FsR_hs = 0.418;               /*! with the relation R_hs/R_ls=1.1. Take the distance between couples of atoms as 1, x = 0.202= 1 - R_hs + R_hs.*/
  system_info.FsR_ls = 0.38;                /*! radio to correspond to low-spin state. */
  system_info.Fsdelta = 0.9;

  system_info.FsTotalAtoms = n * n * n;     /*! total of atoms set, will be FsNi^3, this variable is used like a N (total atoms) */
  system_info.FsNi = n;                     /*! number of atoms per side sqrt(N)*/
  system_info.Fscycles = cycles;            /*! amount of cycles in the Monte Carlos method  */

  system_info.FsScreenrate = 500;
  system_info.FsFilerate = 100;
  system_info.FsSlice = 500;

  system_statistic.sSEnergy_vector.resize(system_info.Fscycles);        /*! sSEnergy_vector, vector that will contain each -Fscycles- energy values to do statistic operations */
  system_statistic.sSLiveAverageEnergy.resize(system_info.Fscycles);
  system_statistic.sSLiveAveE = 0;
  system_statistic.sSFinalAverageEnergy = 0;
  system_statistic.sSMMoment_vector.resize(system_info.Fscycles);       /*! sSMMoment_vector, vector that will contain each -Fscycles- magnetic moment values to do statistic operations */
  system_statistic.sSLiveAverageMMoment.resize(system_info.Fscycles);
  system_statistic.sSLiveAveMM = 0;
  system_statistic.sSFinalAverageMMoment = 0;

  // defining options by default
  system_options.a_temperatureCycle = 0;
  system_options.b_temperatureCycle = 1;
  system_options.totalEnergy_totalMoment = 1;
  system_options.routine = 3;                       /*! call to the short one, without big xyz data.*/
  system_options.metropolisMethod = 2;              /*! call to metropolis with the option one.*/
  system_options.swithEnergy = 3;
  system_options.print_fileout = 1;
  system_options.det_VandL = 1;
  system_options.print_out = 1;
  system_options.liveAverage = 2;
  system_options.xyz_print = 1;
  system_options.temperature_file_print = 1;
  system_options.cycle_file_print = 1;
  system_options.test_function = 2;
  system_options.about_system_inf = 1;
  system_options.VMD_EXPORT_F = 1;
  system_options.open_parammeters_file = 1;

  system_statistic.sSEnergy_vector.resize(system_info.Fscycles);                  //!ver si esta inicialization se hace antes
  std::fill(system_statistic.sSEnergy_vector.begin(), system_statistic.sSEnergy_vector.end(), 0);
  system_statistic.sSLiveAverageEnergy.resize(system_info.Fscycles);              //!
  std::fill(system_statistic.sSLiveAverageEnergy.begin(), system_statistic.sSLiveAverageEnergy.end(), 0);

  system_statistic.sSMMoment_vector.resize(system_info.Fscycles);                 //!ver si esta inicialization se hace antes
  std::fill(system_statistic.sSMMoment_vector.begin(), system_statistic.sSMMoment_vector.end(), 0);
  system_statistic.sSLiveAverageMMoment.resize(system_info.Fscycles);
  std::fill(system_statistic.sSLiveAverageMMoment.begin(), system_statistic.sSLiveAverageMMoment.end(), 0);

  system_statistic.sSEnergy_STD_error_vector.resize(system_info.Fscycles);
  std::fill(system_statistic.sSEnergy_STD_error_vector.begin(), system_statistic.sSEnergy_STD_error_vector.end(), 0);
  system_statistic.sSMoment_STD_error_vector.resize(system_info.Fscycles);
  std::fill(system_statistic.sSMoment_STD_error_vector.begin(), system_statistic.sSMoment_STD_error_vector.end(), 0);

  return 0;
}

int framework::f_define_system_info(framework::Fs_general_info & system_info,
                                    statistic_s & system_statistic,
                                    s_options & system_options)
{
  /*!---------------- initialising system info ----------------*/

  system_info.FsTemperature = 1.0;  /*! system temperature. */
  system_info.FsPressure = 0.01;            /*! system pressure. */
  system_info.FsVolume = 0;                 /*! system volume. */
  system_info.FsL = 0;                      /*! FsL will be modified or calculate. */
  system_info.FsTotalEnergy = 0;            /*! will be modified or calculate. */
  system_info.FsTotalMMoment = 0;           /*! will be modified or calculate. */
  system_info.Fsk_nn = 10;                  /*! Fsk_nn begin in 1 but will be change, but will be a function as soon as possible, is in the order of K_B (e-23). */
  system_info.FsMinT = 1.0;                 /*! min value of temperature cycle. */
  system_info.FsMaxT = 10.0;                  /*! max value of temperature cycle. */
  system_info.FsTsteps = 0.2;
  system_info.FsD = 1.0;                    /*! D represent the value of energy between HS-LS state, is in the order of K_B (e-23). */
  system_info.Fsj_isin = 1.0;              /*! Fsj_isin for this moment represent the value interaction in on a system isin-like mode, is in the order of K_B (e-23). */
  system_info.FsTEnergy_error = 0;          /*! FsTEnergy_error is used as the representation absolute energy error per temperature. */
  system_info.FsTMMoment_error = 0;         /*! FsTMMoment_error is used as the representation absolute magnetic moment error per temperature. */
  system_info.FsR_hs = 0.418;               /*! with the relation R_hs/R_ls=1.1. Take the distance between couples of atoms as 1, x = 0.202= 1 - R_hs + R_hs.*/
  system_info.FsR_ls = 0.38;                /*! radio to correspond to low-spin state. */
  system_info.Fsdelta = 1.0;

  system_info.FsNi = 4;                     /*! number of atoms per side sqrt(N)*/
  system_info.FsTotalAtoms = system_info.FsNi *
                             system_info.FsNi *
                             system_info.FsNi;     /*! total of atoms set, will be FsNi^3, this variable is used like a N (total atoms) */
  system_info.Fscycles = 2000;            /*! amount of cycles in the Monte Carlos method  */

  system_info.FsScreenrate = 500;
  system_info.FsFilerate = 100;
  system_info.FsSlice = 500;

  system_statistic.sSEnergy_vector.resize(system_info.Fscycles);        /*! sSEnergy_vector, vector that will contain each -Fscycles- energy values to do statistic operations */
  system_statistic.sSLiveAverageEnergy.resize(system_info.Fscycles);
  system_statistic.sSLiveAveE = 0;
  system_statistic.sSFinalAverageEnergy = 0;
  system_statistic.sSMMoment_vector.resize(system_info.Fscycles);       /*! sSMMoment_vector, vector that will contain each -Fscycles- magnetic moment values to do statistic operations */
  system_statistic.sSLiveAverageMMoment.resize(system_info.Fscycles);
  system_statistic.sSLiveAveMM = 0;
  system_statistic.sSFinalAverageMMoment = 0;

  // defining options by default
  system_options.a_temperatureCycle = 0;
  system_options.b_temperatureCycle = 1;
  system_options.totalEnergy_totalMoment = 1;
  system_options.routine = 3;                       /*! call to the short one, without big xyz data.*/
  system_options.metropolisMethod = 2;              /*! call to metropolis with the option one.*/
  system_options.swithEnergy = 3;
  system_options.print_fileout = 1;
  system_options.det_VandL = 1;
  system_options.print_out = 1;
  system_options.liveAverage = 2;
  system_options.xyz_print = 1;
  system_options.temperature_file_print = 1;
  system_options.cycle_file_print = 1;
  system_options.test_function = 2;
  system_options.about_system_inf = 1;
  system_options.VMD_EXPORT_F = 1;
  system_options.open_parammeters_file = 1;

  system_statistic.sSEnergy_vector.resize(system_info.Fscycles);                  //!ver si esta inicialization se hace antes
  std::fill(system_statistic.sSEnergy_vector.begin(), system_statistic.sSEnergy_vector.end(), 0);
  system_statistic.sSLiveAverageEnergy.resize(system_info.Fscycles);              //!
  std::fill(system_statistic.sSLiveAverageEnergy.begin(), system_statistic.sSLiveAverageEnergy.end(), 0);

  system_statistic.sSMMoment_vector.resize(system_info.Fscycles);                 //!ver si esta inicialization se hace antes
  std::fill(system_statistic.sSMMoment_vector.begin(), system_statistic.sSMMoment_vector.end(), 0);
  system_statistic.sSLiveAverageMMoment.resize(system_info.Fscycles);
  std::fill(system_statistic.sSLiveAverageMMoment.begin(), system_statistic.sSLiveAverageMMoment.end(), 0);

  system_statistic.sSEnergy_STD_error_vector.resize(system_info.Fscycles);
  std::fill(system_statistic.sSEnergy_STD_error_vector.begin(), system_statistic.sSEnergy_STD_error_vector.end(), 0);
  system_statistic.sSMoment_STD_error_vector.resize(system_info.Fscycles);
  std::fill(system_statistic.sSMoment_STD_error_vector.begin(), system_statistic.sSMoment_STD_error_vector.end(), 0);

  return 0;
}

int framework::create_lattice (std::vector<framework::struct_Sys_properties> & system_in,
	                           std::vector<std::vector<int> > &neighbour_vector,
	                           Fs_general_info & system_info,
	                           const float user_R_hs_d)
{
  std::cout <<"\n----------Creating lattice with "<< system_info.FsTotalAtoms << " atoms.----------"<<std::endl;
  std::cout << "system_info.Fsdelta:" <<system_info.Fsdelta <<std::endl;
  /*!---------------- creating system and neighbour of external variable ----------------*/

  system_in.resize(system_info.FsTotalAtoms);           /*! creating the external vector system.*/
  neighbour_vector.resize(system_info.FsTotalAtoms);    /*! creating the external neighbour vector system.*/

  /*!---------------- initialising system system and neighbour ----------------*/
  //double x_b = double (0.9 - (system_info.FsR_ls + system_info.FsR_ls));
  //f_det_VandL(system_in, system_info, 1);
  int flag = 0;

  for (int z = 0; z < system_info.FsNi; ++z){
    for (int y = 0; y < system_info.FsNi; ++y){
      for (int x = 0; x < system_info.FsNi; ++x){
        system_in[flag].tSpin = 0;
        system_in[flag].sBorder.resize(6); //int_randomNumber_generator(0, 1);
        if (system_in[flag].tSpin == 1) system_in[flag].R = user_R_hs_d; //cambiar esto por una operacion matematica
          else  system_in[flag].R = user_R_hs_d/1.1;
        system_in[flag].spin=1;
        system_in[flag].x =  (2 * x + 1) * (system_info.Fsdelta / 2 + system_info.FsR_hs);
        system_in[flag].y =  (2 * y + 1) * (system_info.Fsdelta / 2 + system_info.FsR_hs);
        system_in[flag].z =  (2 * z + 1) * (system_info.Fsdelta / 2 + system_info.FsR_hs);
        neighbour_vector[flag].resize(6); // 6 vecinos
        flag++;
      }
    }
  }
  std::cout << "system_info.Fsdelta:" <<system_info.Fsdelta <<std::endl;
  //system_info.FsL = f_det_VandL (system_in, system_info, 2);
  system_info.FsL = system_info.Fsdelta * system_info.FsNi;
  std::cout<< "system_info.FsL:"<<system_info.FsL<<std::endl;
  //std::cout << "system_info.FsL\tsystem_info.Fsdelta\tmain_system_info.FsNi\n" <<system_info.FsL <<"\t"<<system_info.Fsdelta <<"\t" <<main_system_info.FsNi<<std::endl;
  system_info.FsVolume = system_info.FsL * system_info.FsL * system_info.FsL;
  std::puts("Done!");
  return 0;
}

int framework::f_temperatureCycle(std::vector<struct_Sys_properties> & system_in,
                                  std::vector<std::vector<int> > & neighbour,
                                  Fs_general_info& system_info,
                                  statistic_s & system_statistic,
                                  s_options & system_options,
                                  const double iniT,
                                  const double finT,
                                  const double Tsteps)
{
  std::cout <<"\n------------------- Temperature cycle (a)-------------------"<<std::endl;
  std::cout <<"\n-------------------------- Options -------------------------"<<std::endl;
  std::cout <<"\n------------------- Cycles:"<<system_options.routine<<" -------------------"<<std::endl;
  std::cout <<"\n------------------- Metropolis:"<<system_options.metropolisMethod<<" -------------------"<<std::endl;
  std::cout <<"\n------------------- Total_Mmoment&Energy:"<<system_options.totalEnergy_totalMoment<<" -------------------"<<std::endl;

  system_info.FsMinT = iniT;
  system_info.FsMaxT = finT;
  int A_TSteps = double(finT - iniT) / Tsteps;

  f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, 0, system_options.totalEnergy_totalMoment);   // OK

  //! to print the xyz coordinates of each atoms before do whatever change. case1 create the head of file and case 2 write the data.
  f_xyz_print(system_in, system_info, system_statistic, 0/*step*/, 1/*!case1*/);
  f_xyz_print(system_in, system_info, system_statistic, 0/*step*/, 2/*!case2*/);

  //! create the head of file_out and the starting point values.
  f_cycle_file_print(system_in, system_info, system_statistic, system_options, 1/*case1*/, 0/*step*/);
  f_cycle_file_print(system_in, system_info, system_statistic, system_options, 2/*case2*/, 0/*step*/);

  //! create the xyz to read in VMD software
  VMD_EXPORT_F(system_in, system_info.FsTotalAtoms, 1/*case*/);

  //! creating the head, this file will contain the data for each temperature.
  f_temperature_file_print(system_in, system_info, system_statistic, 1);

  f_about_system_inf(system_info ,1);

  /**************** TEMPERATURE LOOP ****************/

  for (int i = 0; i < A_TSteps; ++i) {

    system_info.FsTemperature = double (iniT + i * Tsteps);

    std::cout <<"=========== Temperature:"<<system_info.FsTemperature<<" ==========="<<std::endl;
    std::cout <<"=========== Calling the MC cycles =============="<<std::endl;
    std::cout <<"============= Routine option:"<<system_options.routine<<" ================"<<std::endl;

    //! calling routine function with the case 3
    f_routine(system_in, neighbour, system_info, system_statistic, system_options, system_options.routine /*case3*/, system_info.Fscycles);

    f_temperature_file_print(system_in, system_info, system_statistic, 2);
    f_cycle_file_print(system_in, system_info, system_statistic, system_options, 3, system_info.Fscycles - 1);
  }

  return 0;
}

int framework::f_temperatureCycle(std::vector<struct_Sys_properties> & system_in,
                                  std::vector<std::vector<int> > & neighbour,
                                  Fs_general_info& system_info,
                                  statistic_s & system_statistic,
                                  s_options & system_options,
                                  const int opt)
{

  system_statistic.sSEnergy_vector.resize(system_info.Fscycles);                  //!ver si esta inicialization se hace antes
  std::fill(system_statistic.sSEnergy_vector.begin(), system_statistic.sSEnergy_vector.end(), 0);
  system_statistic.sSLiveAverageEnergy.resize(system_info.Fscycles);              //!
  std::fill(system_statistic.sSLiveAverageEnergy.begin(), system_statistic.sSLiveAverageEnergy.end(), 0);

  system_statistic.sSMMoment_vector.resize(system_info.Fscycles);                 //!ver si esta inicialization se hace antes
  std::fill(system_statistic.sSMMoment_vector.begin(), system_statistic.sSMMoment_vector.end(), 0);
  system_statistic.sSLiveAverageMMoment.resize(system_info.Fscycles);
  std::fill(system_statistic.sSLiveAverageMMoment.begin(), system_statistic.sSLiveAverageMMoment.end(), 0);

  system_statistic.sSEnergy_STD_error_vector.resize(system_info.Fscycles);
  std::fill(system_statistic.sSEnergy_STD_error_vector.begin(), system_statistic.sSEnergy_STD_error_vector.end(), 0);
  system_statistic.sSMoment_STD_error_vector.resize(system_info.Fscycles);
  std::fill(system_statistic.sSMoment_STD_error_vector.begin(), system_statistic.sSMoment_STD_error_vector.end(), 0);

  std::cout <<"\n------------------- Temperature cycle -------------------"<<std::endl;
  std::cout <<"\n-----------Temperature range:["<<system_info.FsMinT<<":"<<system_info.FsTsteps<<":"<<system_info.FsMaxT<<"] -----------------------"<<std::endl;
  std::cout <<"\n------------------------ Options: -----------------------"<<std::endl;
  std::cout <<"\n------------------- T_Cycles:"<<system_options.b_temperatureCycle<<" -------------------"<<std::endl;
  std::cout <<"\n------------------- Cycle option:"<<system_options.routine<<" -------------------"<<std::endl;
  std::cout <<"\n------------------- Metropolis option:"<<system_options.metropolisMethod<<" -------------------"<<std::endl;
  std::cout <<"\n------------------- Total_Mmoment&Energy:"<<system_options.totalEnergy_totalMoment<<" -------------------"<<std::endl;

  int A_TSteps = double(system_info.FsMaxT - system_info.FsMinT) / system_info.FsTsteps;


  switch(opt){

    case 1:{ /*!f_temperatureCycle. case.1 General running.*/


      //! restarting vector and statistic variables.
      std::fill(system_statistic.sSEnergy_vector.begin(), system_statistic.sSEnergy_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageEnergy.begin(), system_statistic.sSLiveAverageEnergy.end(), 0);

      std::fill(system_statistic.sSMMoment_vector.begin(), system_statistic.sSMMoment_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageMMoment.begin(), system_statistic.sSLiveAverageMMoment.end(), 0);

      std::fill(system_statistic.sSEnergy_STD_error_vector.begin(), system_statistic.sSEnergy_STD_error_vector.end(), 0);
      std::fill(system_statistic.sSMoment_STD_error_vector.begin(), system_statistic.sSMoment_STD_error_vector.end(), 0);

      f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, 0, system_options.totalEnergy_totalMoment);   // OK

      //! to print the xyz coordinates of each atoms before do whatever change. case1 create the head of file and case 2 write the data.
      f_xyz_print(system_in, system_info, system_statistic, 0/*step*/, 1/*!case1*/);
      f_xyz_print(system_in, system_info, system_statistic, 0/*step*/, 2/*!case2*/);

      //! create the head of file_out and the starting point values.
      f_cycle_file_print(system_in, system_info, system_statistic, system_options, 1/*case1*/, 0/*step*/);
      f_cycle_file_print(system_in, system_info, system_statistic, system_options, 2/*case2*/, 0/*step*/);

      //! create the xyz to read in VMD software
      VMD_EXPORT_F(system_in, system_info.FsTotalAtoms, 1/*case*/);

      //! creating the head, this file will contain the data for each temperature.
      f_temperature_file_print(system_in, system_info, system_statistic, 1);

      f_about_system_inf(system_info ,1);

      /**************** TEMPERATURE LOOP ****************/

      for (int i = 0; i <= A_TSteps; ++i) {

        system_info.FsTemperature = double (system_info.FsMinT + i * system_info.FsTsteps);

        std::cout <<"=========== General case T:"<<system_info.FsTemperature<<"K ==========="<<std::endl;
        std::cout <<"=========== Temperature:"<<system_info.FsTemperature<<" ==========="<<std::endl;
        std::cout <<"=========== Calling the MC cycles =============="<<std::endl;
        std::cout <<"============= Routine option:"<<system_options.routine<<" ================"<<std::endl;
        std::cout <<"============= constant D:"<<system_info.FsD<<" ================"<<std::endl;

        //! calling routine function with the case 3
        f_routine(system_in, neighbour, system_info, system_statistic, system_options, system_options.routine /*case3*/, system_info.Fscycles);

        f_temperature_file_print(system_in, system_info, system_statistic, 2);
        f_cycle_file_print(system_in, system_info, system_statistic, system_options, 3, system_info.Fscycles - 1);
        std::cout<<"STD:"<<system_statistic.sSMoment_STD_error<<std::endl;
      }

      return 0;
      break;
    }

    case 2:{ /*!f_temperatureCycle. case.2 To run the 3D Isin model.*/

      f_remove_olf_files();

      //! restarting vector and statistic variables.
      std::fill(system_statistic.sSEnergy_vector.begin(), system_statistic.sSEnergy_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageEnergy.begin(), system_statistic.sSLiveAverageEnergy.end(), 0);

      std::fill(system_statistic.sSMMoment_vector.begin(), system_statistic.sSMMoment_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageMMoment.begin(), system_statistic.sSLiveAverageMMoment.end(), 0);

      std::fill(system_statistic.sSEnergy_STD_error_vector.begin(), system_statistic.sSEnergy_STD_error_vector.end(), 0);
      std::fill(system_statistic.sSMoment_STD_error_vector.begin(), system_statistic.sSMoment_STD_error_vector.end(), 0);


      //!remove old files.
      std::ofstream isinflux ("isin.dat", std::ios::out|std::ios::app);
       if (!isinflux.is_open()) {
         std::wcerr << "Error abriendo \"isin.dat\" para escribir en el." << std::endl;
         abort ();
      }


      isinflux <<"#step\tTemperature\tAtom\t\tS\tEnergy\tAverEnergy\tSTD_Live_Energy\tMagnetic_Moment\tAveMM\tSTD_live_Moment"<<std::endl;

      for (int i = 0; i <= A_TSteps; ++i) {

        system_info.FsTemperature = double (system_info.FsMinT + i * system_info.FsTsteps);

        std::cout<<std::endl;
        std::cout <<"=========== Temperature(Isin-case): "<<system_info.FsTemperature<<" ==========="<<std::endl;
        std::cout <<"============== Calling the MC cycles =============="<<std::endl;
        std::cout <<"=================== Option(4):"<<system_options.routine<<" ==================="<<std::endl;

        f_routine(system_in, neighbour, system_info, system_statistic, system_options, 4, system_info.Fscycles);

        //f_temperature_file_print(system_in, system_info, system_statistic, 3);
        //f_cycle_file_print(system_in, system_info, system_statistic, system_options, 3, system_info.Fscycles - 1);
      }


      //isinflux <<std::endl;
      isinflux.close();
      return 0;
      break;
    }

    default: std::cout<<"Calling Temperature function with a bad option.\n Type ./vector --help to more information."<<std::endl;
  }


  return 0;
}


/***********************************/
/*                                 */
/*     framework::f_routine ()     */
/*              TIME LOOP          */
/***********************************/

int framework::f_routine(std::vector<struct_Sys_properties> & system_in,
                          const std::vector<std::vector<int> > & neighbour,
                          Fs_general_info & system_info,
                          statistic_s & system_statistic,
                          s_options & system_options,
                          const int opt,
                          const int n_cycles)
{
  switch (opt){

    case 1: { /*!f_routine.opt.1. don't use this fiction yet!*/

      //! restarting vector and statistic variables.
      std::fill(system_statistic.sSEnergy_vector.begin(), system_statistic.sSEnergy_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageEnergy.begin(), system_statistic.sSLiveAverageEnergy.end(), 0);

      std::fill(system_statistic.sSMMoment_vector.begin(), system_statistic.sSMMoment_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageMMoment.begin(), system_statistic.sSLiveAverageMMoment.end(), 0);

      std::fill(system_statistic.sSEnergy_STD_error_vector.begin(), system_statistic.sSEnergy_STD_error_vector.end(), 0);
      std::fill(system_statistic.sSMoment_STD_error_vector.begin(), system_statistic.sSMoment_STD_error_vector.end(), 0);

      std::cout <<"\n------------------- Time cycle (case 1) -------------------"<<std::endl;

      system_statistic.sSEnergy_vector.resize(system_info.Fscycles);
      system_statistic.sSLiveAverageEnergy.resize(system_info.Fscycles);
      system_statistic.sSMMoment_vector.resize(system_info.Fscycles);
      system_statistic.sSLiveAverageMMoment.resize(system_info.Fscycles);

      f_print_fileout(system_in, neighbour, system_info, system_statistic, 2, 0);

      for (int i = 0; i < system_info.Fscycles; ++i){
        std::cout <<"\n=============== "<< i <<" ==============="<<std::endl;
        f_metropolisMethod(system_in, neighbour, system_info, system_statistic, system_options, 2);
        f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, i, 1);
        f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, i, 2);

        if(i%1 == 0){

          VMD_EXPORT_F (system_in, system_info.FsTotalAtoms, 2);
          f_print_fileout(system_in, neighbour, system_info, system_statistic, 3, i);
        }
      }


      break;
    }

    case 2: { /*!f_routine.opt.2-like opt 3. print all files include the big-one.*/

      //! restarting vector and statistic variables.
      std::fill(system_statistic.sSEnergy_vector.begin(), system_statistic.sSEnergy_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageEnergy.begin(), system_statistic.sSLiveAverageEnergy.end(), 0);

      std::fill(system_statistic.sSMMoment_vector.begin(), system_statistic.sSMMoment_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageMMoment.begin(), system_statistic.sSLiveAverageMMoment.end(), 0);

      std::fill(system_statistic.sSEnergy_STD_error_vector.begin(), system_statistic.sSEnergy_STD_error_vector.end(), 0);
      std::fill(system_statistic.sSMoment_STD_error_vector.begin(), system_statistic.sSMoment_STD_error_vector.end(), 0);

      std::cout <<"\n------------------- Time cycle (case 2) -------------------"<<std::endl;


      for (int i = 0; i < system_info.Fscycles; ++i){

        //! calling metropolis function
        f_metropolisMethod(system_in, neighbour, system_info, system_statistic, system_options, system_options.metropolisMethod);                                //!

        //! calculate the energy system and total moment of the system, calling the case 1 and 2 respectively
        f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, i, 1/*case*/);
        f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, i, 2/*case*/);

        //! liveAverage in the case 2 it is use to calculate the average step by step. the same situation with STD error function.
        f_liveAverage(system_statistic, i, 2/*case*/);
        f_STD_desviation (system_info, system_statistic, i/*step*/, 2/*opt*/);

        //! print the principal result of each step.
        f_cycle_file_print(system_in, system_info, system_statistic, system_options, 2, i);
        f_xyz_print(system_in, system_info, system_statistic, i, 2);              /*! big data*/

        if(i % system_info.FsFilerate == 0) VMD_EXPORT_F (system_in, system_info.FsTotalAtoms, 2);
        if(i % system_info.FsScreenrate == 0) std::cout <<"\n=============== "<< i <<" ==============="<<std::endl;
      }
      f_liveAverage(system_statistic, system_info.Fscycles - 1/*step*/, 1/*case*/);
      f_STD_desviation (system_info, system_statistic, system_info.Fscycles - 1/*step*/, 1/*opt*/);

      break;
    }

    case 3: {  /*!f_routine.opt.3. best option.*/

      //! restarting vector and statistic variables.
      std::fill(system_statistic.sSEnergy_vector.begin(), system_statistic.sSEnergy_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageEnergy.begin(), system_statistic.sSLiveAverageEnergy.end(), 0);

      std::fill(system_statistic.sSMMoment_vector.begin(), system_statistic.sSMMoment_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageMMoment.begin(), system_statistic.sSLiveAverageMMoment.end(), 0);

      std::fill(system_statistic.sSEnergy_STD_error_vector.begin(), system_statistic.sSEnergy_STD_error_vector.end(), 0);
      std::fill(system_statistic.sSMoment_STD_error_vector.begin(), system_statistic.sSMoment_STD_error_vector.end(), 0);

      std::cout <<"\n------------------- Time cycle (case 3)-------------------"<<std::endl;

      for (int i = 0; i < system_info.Fscycles; ++i){

        //! calling metropolis function
        f_metropolisMethod(system_in, neighbour, system_info, system_statistic, system_options, system_options.metropolisMethod);                                //!

        //! calculate the energy system and total moment of the system, calling the case 1 and 2 respectively
        f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, i, 1/*case*/);
        f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, i, 2/*case*/);

        //! liveAverage in the case 2 it is use to calculate the average step by step. the same situation with STD error function.
        f_liveAverage(system_statistic, i, 2/*case*/);
        f_STD_desviation (system_info, system_statistic, i/*step*/, 2/*opt*/);

        //! print the principal result of each step.
        f_cycle_file_print(system_in, system_info, system_statistic, system_options, 2, i);

        if(i % system_info.FsFilerate == 0) VMD_EXPORT_F (system_in, system_info.FsTotalAtoms, 2);
        if(i % system_info.FsScreenrate == 0) std::cout <<"\n=============== "<< i <<" ==============="<<std::endl;
      }

      f_liveAverage(system_statistic, system_info.Fscycles - 1/*step*/, 1/*case*/);
      f_STD_desviation(system_info, system_statistic, system_info.Fscycles - 1/*step*/, 1/*opt*/);

      break;
    }

    case 4: {  /*!f_routine.opt.4. running Isin-like model. */

      //! restarting vector and statistic variables.
      std::fill(system_statistic.sSEnergy_vector.begin(), system_statistic.sSEnergy_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageEnergy.begin(), system_statistic.sSLiveAverageEnergy.end(), 0);

      std::fill(system_statistic.sSMMoment_vector.begin(), system_statistic.sSMMoment_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageMMoment.begin(), system_statistic.sSLiveAverageMMoment.end(), 0);

      std::fill(system_statistic.sSEnergy_STD_error_vector.begin(), system_statistic.sSEnergy_STD_error_vector.end(), 0);
      std::fill(system_statistic.sSMoment_STD_error_vector.begin(), system_statistic.sSMoment_STD_error_vector.end(), 0);

      std::cout <<"\n------------- Ising - Time cycle (case 4) ------------"<<std::endl;


      std::ofstream isinflux ("isin.dat", std::ios::out|std::ios::app);
      if (!isinflux.is_open()) {
        std::wcerr << "Error abriendo \"isin.dat\" para escribir en el." << std::endl;
        abort ();
      }

      VMD_EXPORT_F(system_in, system_info.FsTotalAtoms, 2);                             //!
      f_xyz_print(system_in, system_info, system_statistic, 0, 2);

      for (int i = 0; i < system_info.Fscycles; ++i){

        f_metropolisMethod(system_in, neighbour, system_info, system_statistic, system_options, 3);                                //!

        f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, i, 4); //! call the Isin-like model



        f_liveAverage(system_statistic, i, 2);
        f_STD_desviation (system_info, system_statistic, i/*step*/, 2/*opt*/);
        //f_cycle_file_print(system_in, system_info, system_statistic, system_options,2, i);
        isinflux <<i<<"\t"
                 <<system_info.FsTemperature<<"\t"
                 <<"Atom\t"
                 <<"S\t"
                 <<system_statistic.sSEnergy_vector[i]<<"\t"
                 <<system_statistic.sSLiveAverageEnergy[i]<<"\t"
                 <<system_statistic.sSEnergy_STD_error_vector[i]<<"\t"
                 <<system_statistic.sSMMoment_vector[i]<<"\t"
                 <<system_statistic.sSLiveAverageMMoment[i]<<"\t"
                 <<system_statistic.sSMoment_STD_error_vector[i]<<std::endl;

        if(i % system_info.FsFilerate == 0) VMD_EXPORT_F (system_in, system_info.FsTotalAtoms, 2);
        if(i % system_info.FsScreenrate == 0) std::cout <<"\n=============== "<< i <<" ==============="<<std::endl;
      }
      isinflux<<std::endl;
       f_liveAverage(system_statistic, system_info.Fscycles - 1/*step*/, 1/*case*/);
      f_STD_desviation (system_info, system_statistic, system_info.Fscycles - 1/*step*/, 1/*opt*/);

      isinflux.close();



      break;
    }

    case 5: {  /*!f_routine.opt.5. Calculate the energy of spring xyz to print in outut file */

      //! restarting vector and statistic variables.
      std::fill(system_statistic.sSEnergy_vector.begin(), system_statistic.sSEnergy_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageEnergy.begin(), system_statistic.sSLiveAverageEnergy.end(), 0);

      std::fill(system_statistic.sSMMoment_vector.begin(), system_statistic.sSMMoment_vector.end(), 0);
      std::fill(system_statistic.sSLiveAverageMMoment.begin(), system_statistic.sSLiveAverageMMoment.end(), 0);

      std::fill(system_statistic.sSEnergy_STD_error_vector.begin(), system_statistic.sSEnergy_STD_error_vector.end(), 0);
      std::fill(system_statistic.sSMoment_STD_error_vector.begin(), system_statistic.sSMoment_STD_error_vector.end(), 0);

      std::cout <<"\n------------------- Time cycle (case 5)-------------------"<<std::endl;

      for (int i = 0; i < system_info.Fscycles; ++i){

        //! calling metropolis function
        f_metropolisMethod(system_in, neighbour, system_info, system_statistic, system_options, system_options.metropolisMethod);                                //!

        //! calculate the energy system and total moment of the system, calling the case 1 and 2 respectively
        f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, i, system_options.totalEnergy_totalMoment /*case*/);
        f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, i, 2 /*case*/);

        //! liveAverage in the case 2 it is use to calculate the average step by step. the same situation with STD error function.
        f_liveAverage(system_statistic, i, 2/*case*/);
        f_STD_desviation (system_info, system_statistic, i/*step*/, 2/*opt*/);

        //! print the principal result of each step.
        f_cycle_file_print(system_in, system_info, system_statistic, system_options, 2, i);

        if(i % system_info.FsFilerate == 0) VMD_EXPORT_F (system_in, system_info.FsTotalAtoms, 2);
        if(i % system_info.FsScreenrate == 0) std::cout <<"=============== "<< i <<" ==============="<<std::endl;
      }

      f_liveAverage(system_statistic, system_info.Fscycles - 1/*step*/, 1/*case*/);
      f_STD_desviation(system_info, system_statistic, system_info.Fscycles - 1/*step*/, 1/*opt*/);

      break;
    }

    default: printf("framework::f_routine. Choice other than 1, 2, 3 or 4.\n");
    break;
  }
  return 0;
}


/******************************************************/
/*                                                    */
/*      framework::f_totalEnergy_totalMoment ()       */
/* TO CALCULATE THE TOTAL MOMENT AND THE TOTAL ENERGY */
/*                                                    */
/******************************************************/

double framework::f_totalEnergy_totalMoment (std::vector<struct_Sys_properties> & system_in,
	                                           const std::vector<std::vector<int> > &neighbour,
	                                           Fs_general_info & system_info,
                                             statistic_s & system_statistic,
                                             const int index,
	                                           const int opt)
{
  //std::cout<<"---------------------------- "<<opt<<" ----------------------------"<<std::endl;
  double DOUBLE_E_OR_MM,
         H_o = 0,
         H_nn = 0,
         H_mag = 0,
         dx = 0,
         dy = 0,
         dz = 0,
         r = 0,
         Total_MMoment= 0;

  switch (opt){

    case 1: { /*!opt.1 (Total_Energy)*/

      system_info.FsTotalEnergy = 0;

      /*! estimating the H_o */

      for (int i = 0; i < system_info.FsTotalAtoms; ++i) H_o += system_in[i].tSpin;
      H_o *= system_info.FsD;



      /*! estimating the H_nn */
      for (int e = 0; e < system_info.FsTotalAtoms; ++e){
        for (int i = 0; i < 6; ++i){
          r = f_tow_p_distance(system_in, system_info, neighbour, e, i);
          H_nn += (r - (system_in[e].R + system_in[neighbour[e][i]].R)) *
                  (r - (system_in[e].R + system_in[neighbour[e][i]].R));
        }
      }
      H_nn *= (system_info.Fsk_nn / 2);



      /*! estimating the H_mag */
      for (int e = 0; e < system_info.FsTotalAtoms; ++e){
        for (int i = 0; i < 6; ++i) H_mag += system_in[e].tSpin * system_in[ neighbour[e][i]].tSpin;
      }
      H_mag *= (double) - system_info.Fsj_isin;

      system_statistic.sSEnergy_vector[index] = H_o + H_nn + H_mag;

      return system_statistic.sSEnergy_vector[index];

      break;
    }

    case 2: { /*!opt.2 (Total_Moment)*/
      double Total_MMoment = 0;
      system_info.FsTotalMMoment = 0;

      for (int i = 0; i < system_info.FsTotalAtoms; ++i) Total_MMoment += system_in[i].tSpin;

      system_statistic.sSMMoment_vector[index] = Total_MMoment/system_info.FsTotalAtoms;// / system_info.FsTotalAtoms;

      /*
      std::cout <<"----------------- lattice-Metro-despues:"<<system_in.size() <<" --------------------"<<std::endl;
      int bandera = 0;
      for (int i = 0; i < system_info.FsTotalAtoms/system_info.FsNi; ++i){
        for (int e = 0; e < system_info.FsNi; ++e){
          std::cout<< system_in[bandera].tSpin <<"\t" ;
          bandera ++;
        }
        std::cout<< std::endl;
      }
      */

      return system_statistic.sSMMoment_vector[index];
      break;
    }

    case 3: { /*!opt.3 Energy of the sprint:(xyz)*/
      system_info.FsTotalEnergy = 0;
      /*! estimating the H_o */
      for (int i = 0; i < system_info.FsTotalAtoms; ++i) H_o += system_in[i].tSpin;
      H_o *= system_info.FsD;

      /*! estimating the H_nn */
      for (int e = 0; e < system_info.FsTotalAtoms; ++e){
        for (int i = 0; i < 6; ++i){
          r = f_tow_p_distance(system_in, system_info, neighbour, e, i);
          H_nn += (r - (system_in[e].R + system_in[neighbour[e][i]].R)) *
                  (r - (system_in[e].R + system_in[neighbour[e][i]].R));
        }
      }
      H_nn *= (system_info.Fsk_nn / 2);

      system_statistic.sSEnergy_vector[index] = H_o + H_nn;

      return system_statistic.sSEnergy_vector[index];
      break;
    }

    case 4: { /*!opt.4(Isin like model)*/
      //std::cout <<"\n------------------- Determining the Isin-like model  -------------------"<<std::endl;


      /*!----------- estimating the H_nn -------------*/
      for (int e = 0; e < system_info.FsTotalAtoms; ++e){
        for (int i = 0; i < 6; ++i){
          H_mag += system_in[e].spin * system_in[ neighbour[e][i]].spin;
        }
        Total_MMoment += system_in[e].spin;   /*! summation over all atoms to estimate the total magnetic moments*/
      }

      H_mag *= system_info.Fsj_isin;
      system_statistic.sSEnergy_vector[index] = - H_mag;
      system_statistic.sSMMoment_vector[index] = Total_MMoment/system_info.FsTotalAtoms ;

      return system_statistic.sSEnergy_vector[index];

      break;
    }

    case 5: { /*!opt.5 Energy of spin stabilization*/

      H_mag = 0;
      Total_MMoment = 0;
      for (int e = 0; e < system_info.FsTotalAtoms; ++e){
        for (int i = 0; i < 6; ++i) H_mag += system_in[e].tSpin * system_in[ neighbour[e][i]].tSpin;
      }
      //std::cout<<"suma:"<<H_mag<<std::endl;
      H_mag *= system_info.Fsj_isin;
      system_statistic.sSEnergy_vector[index] = double(-H_mag);

      return system_statistic.sSEnergy_vector[index];
      break;
    }

    default: printf("f_totalEnergy_totalMoment. Choice other than 1 (Total Energy (Momment and XYZ)), 2 (Total Moment), 3(Total Energy xyz) and 4(Isin-like model)\n");
              break;
   }

  return DOUBLE_E_OR_MM;
}

/*****************************************************************************/
/* f_metropolisMethod(): public function to be use externally and internally */
/*                                                                           */
/*****************************************************************************/

int framework::f_metropolisMethod(std::vector<framework::struct_Sys_properties> & system_in,
	                                const std::vector<std::vector<int> > & neighbour,
	                                Fs_general_info & system_info,
                                  statistic_s & system_statistic,
                                  s_options & system_options,
                                  const int opt)
{
  //std::cout <<"\n------------------- metropolis Method -------------------"<<std::endl;
  const double kboltzmann = 1.3807e-23;
  std::vector<int> rPositions;
  rPositions.resize(system_in.size());                                   //var to storage the random position generator to do the metropolis cycle.
  framework::struct_Sys_properties tmp_system;


  //! generate the vector with random position to select the atom to do the calculus.
  for (int i = 0; i < (system_in.size()); ++i) rPositions[i] = i;
  randomSamplePosition(rPositions,system_info.FsNi);                     //call the random position generator

  double deltaEnergySwitch = 0.0;

  switch(opt){

    case 1:{/*! f_metropolisMethod: case 1. running independently the stabilization of spin/xyz position.(1/0)*/

      for (int i = 0; i < (system_in.size()); ++i){

        //! call f_switchEnergy with case 1, to calculate only the spin energy.
        deltaEnergySwitch = f_switchEnergy(system_in,neighbour, tmp_system, system_info, rPositions[i], 1/*case*/);

        // generate the random number and calculate the probability using a Boltzmann distribution.
        double deltaW = deltaEnergySwitch;
	//+ system_info.FsPressure * (system_info.FsVolume - system_info.FsVolume ) - system_info.FsTotalAtoms * system_info.FsTemperature * log((system_info.FsVolume) / (system_info.FsVolume));
        double chi_prob_gen = real_randomNumber_generator(0, 1);
        double proB_cal = exp(-deltaW/ ( system_info.FsTemperature ) );

        if (deltaW <= 0.0)
	{
          system_in[rPositions[i]].R = tmp_system.R;
          system_in[rPositions[i]].tSpin = tmp_system.tSpin;
          proB_cal = 1.0;
        }
        else if (proB_cal > chi_prob_gen) {
          system_in[rPositions[i]].R = tmp_system.R;
          system_in[rPositions[i]].tSpin = tmp_system.tSpin;
        }

        /*std::cout <<"----------------- lattice-Metro-despues:"<<system_in.size() <<" --------------------"<<std::endl;
         bandera = 0;
         for (int i = 0; i < system_info.FsTotalAtoms/system_info.FsNi; ++i){
          for (int e = 0; e < system_info.FsNi; ++e){
            std::cout<< system_in[bandera].tSpin <<"\t" ;
            bandera ++;
          }
          std::cout<< std::endl;
         }
        */

        //! call f_switchEnergy with case 2, to calculate only the xyz stability.
        tmp_system.x = 0;
        tmp_system.y = 0;
        tmp_system.z = 0;

        //if (rPositions[i] == 0)  deltaEnergySwitch = f_switchEnergy(system_in,neighbour, tmp_system, system_info, rPositions[i], 7);
        //else {
        deltaEnergySwitch = f_switchEnergy(system_in,neighbour, tmp_system, system_info, rPositions[i], 2);
        //}

        deltaW = deltaEnergySwitch;// + system_info.FsPressure * (system_info.FsVolume - system_info.FsVolume ) - system_info.FsTotalAtoms * system_info.FsTemperature * log((system_info.FsVolume) / (system_info.FsVolume)); //this is the moss important part of the function, here it is define the energy criteria.....CARE

        chi_prob_gen = real_randomNumber_generator(0, 1);
        proB_cal = exp(-deltaW/ ( system_info.FsTemperature ) );
        //std::cout << "deltaW\tproB_cal\n"<<deltaW<<"\t"<<proB_cal<<std::endl;
        if (deltaW <= 0.0) {
          system_in[rPositions[i]].x = tmp_system.x;
          system_in[rPositions[i]].y = tmp_system.y;
          system_in[rPositions[i]].z = tmp_system.z;
          system_in[rPositions[i]].tSpin = tmp_system.tSpin;
          system_in[rPositions[i]].R = tmp_system.R;
          proB_cal=1.0;
        }
        else if (proB_cal > chi_prob_gen) {
          system_in[rPositions[i]].x = tmp_system.x;
          system_in[rPositions[i]].y = tmp_system.y;
          system_in[rPositions[i]].z = tmp_system.z;
          system_in[rPositions[i]].tSpin = tmp_system.tSpin;
          system_in[rPositions[i]].R = tmp_system.R;
        }
      }
      break;
    }

    case 2:{ /*!f_metropolisMethod: case 2. running and do both changes, spin and xyz position at the same time.(1/0)*/

      for (int i = 0; i < (system_in.size()); ++i){

        /*std::cout <<"----------------- lattice-Metro-antes:"<<system_in.size() <<" --------------------"<<std::endl;
          int bandera = 0;
          for (int i = 0; i < system_info.FsTotalAtoms/system_info.FsNi; ++i){
            for (int e = 0; e < system_info.FsNi; ++e){
            std::cout<< system_in[bandera].tSpin <<"\t" ;
            bandera ++;
          }
          std::cout<< std::endl;
         }
        */

      tmp_system.x = 0;
      tmp_system.y = 0;
      tmp_system.z = 0;

      //! call f_switchEnergy with case 3, to calculate the energy to change the spin and the xyz position.
      deltaEnergySwitch = f_switchEnergy(system_in,neighbour, tmp_system, system_info, rPositions[i], 3/*case*/);

      // generate the random number and calculate the probability using a Boltzmann distribution.
      double deltaW = deltaEnergySwitch;// + system_info.FsPressure * (system_info.FsVolume - system_info.FsVolume ) - system_info.FsTotalAtoms * system_info.FsTemperature * log((system_info.FsVolume) / (system_info.FsVolume)); //this is the moss important part of the function, here it is define the energy criteria.....CARE
      double chi_prob_gen = real_randomNumber_generator(0, 1);
      double proB_cal = exp(-deltaW/ ( system_info.FsTemperature ) );

      if (deltaW <= 0.0)
      {
        system_in[rPositions[i]].x = tmp_system.x;
        system_in[rPositions[i]].y = tmp_system.y;
        system_in[rPositions[i]].z = tmp_system.z;
        system_in[rPositions[i]].R = tmp_system.R;
        system_in[rPositions[i]].tSpin = tmp_system.tSpin;
        //std::cout<<"deltaW:"<<1<<std::endl;
        proB_cal=1.0;
      }
      else if (proB_cal > chi_prob_gen) {
        system_in[rPositions[i]].x = tmp_system.x;
        system_in[rPositions[i]].y = tmp_system.y;
        system_in[rPositions[i]].z = tmp_system.z;
        system_in[rPositions[i]].R = tmp_system.R;
        system_in[rPositions[i]].tSpin = tmp_system.tSpin;
        //std::cout<<"proB_cal > chi_prob_gen:"<<1<<std::endl;
      }

      /*std::cout <<"----------------- lattice-Metro-despues:"<<system_in.size() <<" --------------------"<<std::endl;
        bandera = 0;
        for (int i = 0; i < system_info.FsTotalAtoms/system_info.FsNi; ++i){
          for (int e = 0; e < system_info.FsNi; ++e){
            std::cout<< system_in[bandera].tSpin <<"\t" ;
            bandera ++;
          }
          std::cout<< std::endl;
        }*/
      }
      break;
    }

    case 3:{/*! f_metropolisMethod: case 3. running the Isin model 3D(-1/1)*/

      for (int i = 0; i < (system_in.size()); ++i){

        //! Isin model switch case 4.
        deltaEnergySwitch = f_switchEnergy(system_in,neighbour, tmp_system, system_info, rPositions[i], 4);

        //!
        double chi_prob_gen = real_randomNumber_generator(0, 1);
        double proB_cal = exp(- deltaEnergySwitch/ (system_info.FsTemperature));

        if (deltaEnergySwitch <= 0.0) {
          system_in[rPositions[i]].spin = tmp_system.spin;
          proB_cal=1.0;
        }
        else if (proB_cal > chi_prob_gen) system_in[rPositions[i]].spin = tmp_system.spin;
      }

      break;
    }

    case 4:{/*! f_metropolisMethod: case 4. running only taking in account the spin (1/0).*/

      for (int i = 0; i < (system_in.size()); ++i){

        //for (int e = 0; e <  (system_in.size()); ++e)
        //{
        //  std::cout<<system_in[e].tSpin<<"\t";
        //}
        //std::cout<<std::endl;



        deltaEnergySwitch = f_switchEnergy(system_in,neighbour, tmp_system, system_info, rPositions[i], 1);

        double deltaW = deltaEnergySwitch;
	// + system_info.FsPressure * (system_info.FsVolume - system_info.FsVolume ) - system_info.FsTotalAtoms * system_info.FsTemperature * log((system_info.FsVolume) / (system_info.FsVolume));
        double chi_prob_gen = real_randomNumber_generator(0, 1);
        double proB_cal = exp(-deltaEnergySwitch/ ( system_info.FsTemperature ) );

        if (deltaEnergySwitch <= 0.0) {
          system_in[rPositions[i]].tSpin = tmp_system.tSpin;
         // proB_cal = 1.0;
        }
        else if (proB_cal > chi_prob_gen) system_in[rPositions[i]].tSpin = tmp_system.tSpin;

      }
      break;
    }

    case 5:{/*! f_metropolisMethod: case 5. running only taking in account the xyz stabilization.*/

      for (int i = 0; i < (system_in.size()); ++i){

        tmp_system.x = 0;
        tmp_system.y = 0;
        tmp_system.z = 0;

        //if (rPositions[i] == 0)  deltaEnergySwitch = f_switchEnergy(system_in,neighbour, tmp_system, system_info, rPositions[i], 7);
        //else {
        deltaEnergySwitch = f_switchEnergy(system_in,neighbour, tmp_system, system_info, rPositions[i], 2);
        //}

        double deltaW = deltaEnergySwitch;// + system_info.FsPressure * (system_info.FsVolume - system_info.FsVolume ) - system_info.FsTotalAtoms * system_info.FsTemperature * log((system_info.FsVolume) / (system_info.FsVolume)); //this is the moss important part of the function, here it is define the energy criteria.....CARE

        double chi_prob_gen = real_randomNumber_generator(0, 1);
        double proB_cal = exp(-deltaW/ ( system_info.FsTemperature ) );
        //std::cout << "deltaW\tproB_cal\n"<<deltaW<<"\t"<<proB_cal<<std::endl;
        if (deltaW <= 0.0) {
          system_in[rPositions[i]].x = tmp_system.x;
          system_in[rPositions[i]].y = tmp_system.y;
          system_in[rPositions[i]].z = tmp_system.z;
          system_in[rPositions[i]].tSpin = tmp_system.tSpin;
          system_in[rPositions[i]].R = tmp_system.R;
          proB_cal=1.0;
        }
        else if (proB_cal > chi_prob_gen) {
          system_in[rPositions[i]].x = tmp_system.x;
          system_in[rPositions[i]].y = tmp_system.y;
          system_in[rPositions[i]].z = tmp_system.z;
          system_in[rPositions[i]].tSpin = tmp_system.tSpin;
          system_in[rPositions[i]].R = tmp_system.R;
        }
      }
      break;
    }
    default: std::cout<<"f_metropolisMethod switch"<<std::endl;
  }

  double tmp_L = system_info.FsL;
  double tmp_V = system_info.FsVolume;

  tmp_L = system_info.FsL + float(0.005 * real_randomNumber_generator(-1 ,1) * system_info.FsL);
  tmp_V = tmp_L * tmp_L * tmp_L;

  double deltaW = system_info.FsPressure * (tmp_V - system_info.FsVolume) - system_info.FsTotalAtoms * system_info.FsTemperature * log((tmp_V) / (system_info.FsVolume)); //this is the moss important part of the function, here it is define the energy criteria.....CARE

  double chi_prob_gen = real_randomNumber_generator(0, 1);
  double proB_cal = exp(-deltaW/ ( system_info.FsTemperature ) );

  if (tmp_L <= system_info.FsL && tmp_L >= ((system_info.FsNi * 2 * system_info.FsR_ls))) {
    //std::cout<<"acept new L:"<<tmp_L<<"\tNi:"<<system_info.FsNi<<"\tFsls:"<< system_info.FsR_ls<<std::endl;
    for (size_t i = 0; i < system_info.FsTotalAtoms ; i++) {
      system_in[i].x *= tmp_L / system_info.FsL;
      system_in[i].y *= tmp_L / system_info.FsL;
      system_in[i].z *= tmp_L / system_info.FsL;
    }
    system_info.FsL = tmp_L;
    system_info.FsVolume = system_info.FsL * system_info.FsL * system_info.FsL;
  }
  else if (proB_cal > chi_prob_gen){
    //std::cout<<"acept new L:"<<tmp_L<<"\tNi:"<<system_info.FsNi<<"\tFsls:"<< system_info.FsR_ls<<std::endl;
    for (size_t i = 0; i < system_info.FsTotalAtoms ; i++) {
      system_in[i].x *= tmp_L / system_info.FsL;
      system_in[i].y *= tmp_L / system_info.FsL;
      system_in[i].z *= tmp_L / system_info.FsL;
    }
    system_info.FsL = tmp_L;
    system_info.FsVolume = system_info.FsL * system_info.FsL * system_info.FsL;
  }

  return 0;
}

/**************************************************************/
/* swithEnergy(): public function to return the energy switch */
/*                                                            */
/**************************************************************/

double framework::f_switchEnergy(const std::vector<framework::struct_Sys_properties> & system_in,
	                            const std::vector<std::vector<int> > & neighbour,
	                            framework::struct_Sys_properties & temp_var,
	                            Fs_general_info & system_info,
	                            const int index,
                              const int opt)
{
  double dx, dy, dz, sx, sy, sz, r,
         E_initial = 0,
         E_final = 0,
         H_mag_ini = 0,
         H_mag_fin = 0;

  switch (opt){
    case 1:{ /*!case 1: only for compute the switch energy spin(0/1). -D\Sigma s_i */

      system_info.FsD = ((-system_info.FsTemperature / 5) + 1);

      //! calculate the energy of the actual state.
      H_mag_ini = (double)(-system_info.FsD * system_in[index].tSpin);

      //!changing the configuration and putting in a temporal storage.
      temp_var.tSpin = (1 - system_in[index].tSpin) * 1;
      temp_var.R = (1 - system_in[index].tSpin) * system_info.FsR_hs + system_in[index].tSpin * system_info.FsR_ls;

      //! recalculate the energy of the new state.
      H_mag_fin = (double)(system_info.FsD * temp_var.tSpin);
      //std::cout<<"delta energy:"<< H_mag_fin - H_mag_ini<< std::endl;
      return double(H_mag_fin - H_mag_ini);
    }

    case 2:{/*!case 2: only for compute the switch energy xyz changes k_nn\Sigma (r_i - (R_i + R_j))^2.*/

      //! calculate the energy of the actual state.
      for (int i = 0; i < 6; ++i){
        r = f_tow_p_distance(system_in, system_info, neighbour, index, i);
        E_initial += (r - (system_in[index].R + system_in[neighbour[index][i]].R)) *
                     (r - (system_in[index].R + system_in[neighbour[index][i]].R));
      }
      E_initial *= double(system_info.Fsk_nn / 2);

      //!changing the configuration and putting in a temporal storage.
      temp_var.x = system_in[index].x + float(0.005 * real_randomNumber_generator(-1 ,1) * system_info.FsL);
      temp_var.y = system_in[index].y + float(0.005 * real_randomNumber_generator(-1 ,1) * system_info.FsL);
      temp_var.z = system_in[index].z + float(0.005 * real_randomNumber_generator(-1 ,1) * system_info.FsL);

      temp_var.tSpin = (1 - system_in[index].tSpin) * 1;
      temp_var.R = (1 - system_in[index].tSpin) * system_info.FsR_hs + system_in[index].tSpin * system_info.FsR_ls;

      //! recalculate the energy of the new state.
      for (int i = 0; i < 6; ++i){
        r = f_tow_p_distance(system_in, temp_var , system_info, neighbour, index, i);
        if (r >= 2 * system_info.FsR_hs)
          E_final += (r - (temp_var.R + system_in[neighbour[index][i]].R)) *
                     (r - (temp_var.R + system_in[neighbour[index][i]].R));
        else {
          temp_var.x = system_in[index].x;
          temp_var.y = system_in[index].y;
          temp_var.z = system_in[index].z;
          temp_var.tSpin = system_in[index].tSpin;
          temp_var.R = system_in[index].R;
          return 0;
        }
      }
      E_final *= double(system_info.Fsk_nn / 2);

      return double(E_final - E_initial);
    }

    case 3:{/*!case 3: compute the switch energy doing to both xyz and spin(0/1) changes at the same time.*/

      //! calculate the energy of the actual state.
      for (int i = 0; i < 6; ++i){
        r = f_tow_p_distance(system_in, system_info, neighbour, index, i);
        E_initial += (r - (system_in[index].R + system_in[neighbour[index][i]].R)) *
                     (r - (system_in[index].R + system_in[neighbour[index][i]].R));
      }

      H_mag_ini = double (system_info.FsD * system_in[index].tSpin);

      //!changing the configuration and putting in a temporal storage.
      E_initial *= double(system_info.Fsk_nn / 2);
      H_mag_ini *= (double) - system_info.Fsj_isin;

      sx = system_in[index].x;
      sy = system_in[index].y;
      sz = system_in[index].z;

      temp_var.x = sx + float (0.005 * real_randomNumber_generator(-1 ,1));
      temp_var.y = sy + float (0.005 * real_randomNumber_generator(-1 ,1));
      temp_var.z = sz + float (0.005 * real_randomNumber_generator(-1 ,1));

      temp_var.tSpin = (1 - system_in[index].tSpin) * 1 + system_in[index].tSpin * 0;
      temp_var.R = (1 - system_in[index].tSpin) * system_info.FsR_hs + system_in[index].tSpin * system_info.FsR_ls;

      //! recalculate the energy of the new state.
      for (int i = 0; i < 6; ++i){
        r = f_tow_p_distance(system_in, temp_var , system_info, neighbour, index, i);
        E_final += (r - (temp_var.R + system_in[ neighbour[index][i] ].R)) *
                   (r - (temp_var.R + system_in[ neighbour[index][i] ].R));
        H_mag_fin += temp_var.tSpin + system_in[ neighbour[index][i] ].tSpin;
      }
      E_final *= double(system_info.Fsk_nn / 2);
      H_mag_fin *= (double) -system_info.Fsj_isin;

      return double((E_final + H_mag_fin) - (E_initial + H_mag_ini));
      break;
    }

    case 4:{/*!case 4: compute the switch energy for Isin model (-1/1).*/

      //! calculate the energy of the actual state.
      for (int i = 0; i < 6; ++i) H_mag_ini += system_in[index].spin * system_in[neighbour[index][i]].spin;
      H_mag_ini *= (double) - system_info.Fsj_isin;

      //!changing the configuration and putting in a temporal storage.
      temp_var.spin = system_in[index].spin * -1;

      //! recalculate the energy of the new state.
      for (int i = 0; i < 6; ++i) H_mag_fin += temp_var.spin * system_in[ neighbour[index][i] ].spin;
      H_mag_fin *= (double) - system_info.Fsj_isin;

      return double(H_mag_fin - H_mag_ini);
      break;
    }

    case 5:{

      double H_nnn_ini = 0, H_nnn_fin = 0;
      //! calculate the energy of the actual state.
      for (int i = 0; i < 6; ++i) H_mag_ini += system_in[index].tSpin * system_in[neighbour[index][i]].tSpin;
      H_mag_ini *= (double) - system_info.Fsj_isin;

      for (int i = 0; i < 6; ++i) H_nnn_ini += system_in[index].tSpin * system_in[neighbour[neighbour[index][i]][i]].tSpin;
      H_nnn_ini *= (double) - system_info.Fsj_isin;


      //!changing the configuration and putting in a temporal storage.
      temp_var.tSpin = (1 - system_in[index].tSpin) * 1;
      temp_var.R = (1 - system_in[index].tSpin) * system_info.FsR_hs + system_in[index].tSpin * system_info.FsR_ls;

      //! recalculate the energy of the new state.
      for (int i = 0; i < 6; ++i) H_mag_fin += temp_var.tSpin * system_in[ neighbour[index][i] ].tSpin;
      H_mag_fin *= (double) - system_info.Fsj_isin;

      for (int i = 0; i < 6; ++i) H_nnn_fin += temp_var.tSpin * system_in[neighbour[neighbour[index][i]][i]].tSpin;
      H_nnn_fin *= (double) - system_info.Fsj_isin;

      return double(H_mag_fin + H_nnn_fin ) - (H_mag_ini +  H_nnn_fin);
    }

    case 6:{ /*!case 6: only for compute the switch energy spin(0/1). -D\Sigma s_i*/

     //! calculate the energy of the actual state.
      for (int i = 0; i < 6; ++i) H_mag_ini += system_in[index].tSpin * system_in[neighbour[index][i]].tSpin;
      H_mag_ini *= (double) - system_info.FsD ;

      //!changing the configuration and putting in a temporal storage.
      temp_var.tSpin = (1 - system_in[index].tSpin) * 1;
      temp_var.R = (1 - system_in[index].tSpin) * system_info.FsR_hs + system_in[index].tSpin * system_info.FsR_ls;

      //! recalculate the energy of the new state.
      for (int i = 0; i < 6; ++i) H_mag_fin += temp_var.tSpin * system_in[ neighbour[index][i] ].tSpin;
      H_mag_fin *= (double) - system_info.FsD;

      return double(H_mag_fin - H_mag_ini);
    }

    case 7:{/*!case 7: only for compute the switch energy xyz special case of position 0 changes k_nn\Sigma (r_i - (R_i + R_j))^2.*/

      temp_var.x = system_in[index].x;
      temp_var.y = system_in[index].y;
      temp_var.z = system_in[index].z;

      //! calculate the energy of the actual state.
      for (int i = 0; i < 6; ++i){
        r = f_tow_p_distance(system_in, system_info, neighbour, index, i);
        E_initial += (r - (system_in[index].R + system_in[neighbour[index][i]].R)) *
                     (r - (system_in[index].R + system_in[neighbour[index][i]].R));
      }
      E_initial *= double(system_info.Fsk_nn / 2);



      temp_var.tSpin = (1 - system_in[index].tSpin) * 1;
      temp_var.R = (1 - system_in[index].tSpin) * system_info.FsR_hs + system_in[index].tSpin * system_info.FsR_ls;

      //! recalculate the energy of the new state.
      for (int i = 0; i < 6; ++i){
        //r = f_tow_p_distance(system_in, temp_var , system_info, neighbour, index, i);
        E_final += (r - (temp_var.R + system_in[neighbour[index][i]].R)) *
                   (r - (temp_var.R + system_in[neighbour[index][i]].R));
      }
      E_final *= double(system_info.Fsk_nn / 2);

      return double(E_final - E_initial);
    }

    default: std::cout<<"f_switchEnergy the options available are 1, 2 and 3."<<std::endl;
  }
  //std::cout<<"E_final:"<<E_final<<"\t"<<"E_initial:"<< E_initial<<"\tr:"<<r <<std::endl;
  return (E_final + H_mag_fin) - (E_initial + H_mag_ini) ;
}

/*******************************************************************/
/* framework::f_det_VandL(): public function to return the volume */
/*        and write in the Fs_general_info structure the L         */
/*******************************************************************/

double framework::f_det_VandL (const std::vector<framework::struct_Sys_properties> & system_in,
	                             Fs_general_info & system_info,
	                             const int opt)
{
  switch (opt){
    case 1:{ /*! Determining de value of L and Volume (frist time.opt_1) */
      system_info.FsL = 0;
      system_info.FsVolume = 0;
      //std::cout<<"Determining de value of L and Volume (frist time.opt_1)"<<std::endl;
      system_info.FsL = 0.202 + 2 * system_info.FsR_hs * system_info.FsNi + (system_info.FsNi - 1) * 0.202; //0.202 es la separacion entre la superficie de cada esfera

      system_info.FsVolume  = system_info.FsL * system_info.FsL * system_info.FsL;

      break;
    }

    case 2: { /*!f_det_VandL case2 Determining de value of L  (recalculation.opt_2)*/
      //std::cout<<"Determining de value of L  (recalculation.opt_2)"<<std::endl;
      double dx = 0, dy = 0, dz = 0,
             rx = 0, ry = 0, rz = 0,
             distance = 0;


      dx = (system_in[0].x - system_in[system_info.FsNi - 1].x) * (system_in[0].x - system_in[system_info.FsNi - 1].x);
      dy = (system_in[0].y - system_in[system_info.FsNi - 1].y) * (system_in[0].y - system_in[system_info.FsNi - 1].y);
      dz = (system_in[0].z - system_in[system_info.FsNi - 1].z) * (system_in[0].z - system_in[system_info.FsNi - 1].z);
      rx = sqrt(dx + dy + dz);

      dx = (system_in[0].x - system_in[(system_info.FsNi * system_info.FsNi) - system_info.FsNi].x) * (system_in[0].x - system_in[(system_info.FsNi * system_info.FsNi) - system_info.FsNi].x);
      dy = (system_in[0].y - system_in[(system_info.FsNi * system_info.FsNi) - system_info.FsNi].y) * (system_in[0].y - system_in[(system_info.FsNi * system_info.FsNi) - system_info.FsNi].y);
      dz = (system_in[0].z - system_in[(system_info.FsNi * system_info.FsNi) - system_info.FsNi].z) * (system_in[0].z - system_in[(system_info.FsNi * system_info.FsNi) - system_info.FsNi].z);
      ry = sqrt(dx + dy + dz);

      dx = (system_in[0].x - system_in[system_info.FsTotalAtoms - (system_info.FsNi * system_info.FsNi)].x) * (system_in[0].x - system_in[system_info.FsTotalAtoms - (system_info.FsNi * system_info.FsNi)].x);
      dy = (system_in[0].y - system_in[system_info.FsTotalAtoms - (system_info.FsNi * system_info.FsNi)].y) * (system_in[0].y - system_in[system_info.FsTotalAtoms - (system_info.FsNi * system_info.FsNi)].y);
      dz = (system_in[0].z - system_in[system_info.FsTotalAtoms - (system_info.FsNi * system_info.FsNi)].z) * (system_in[0].z - system_in[system_info.FsTotalAtoms - (system_info.FsNi * system_info.FsNi)].z);

      rz = sqrt(dx + dy + dz);


      return (double) ((rx + ry + rz)/3);

      break;
    }

    default: printf("framework::f_routine. Choice other than 1, 2 and 3\n");
              break;
   }
  return system_info.FsVolume;
}

int framework::f_print_out(const std::vector<struct_Sys_properties> & system_in,
	                         const std::vector<std::vector<int> > & neighbour,
	                         const framework::Fs_general_info & system_info,
	                         const int opt)
{
  switch (opt){
    case 1: { /*!f_print_out.opt.1 (use only to create the head and remove the old files. )*/

      std::cout <<"\n---------- Removing old files to create a new one!. ----------"<<std::endl;



      std::ofstream simple_result_vs_T ("result_vs_T.dat", std::ios::out|std::ios::app);
      if (!simple_result_vs_T.is_open()) {
        std::wcerr << "Error abriendo \"result_vs_T.dat\" para escribir en el." << std::endl;
        abort ();
      }

      simple_result_vs_T <<"#No_Atoms:                         "<<system_info.FsTotalAtoms
                         <<"\n#Temperature range [min:max]:["<<system_info.FsMinT<<":"<<system_info.FsMaxT<<"]"
                         <<"\n#Pressure:                       "<<system_info.FsPressure
                         <<"\n#Volume:                         "<<system_info.FsVolume
                         <<"\n#K_nn:                           "<<system_info.Fsk_nn
                         <<"\n#Total_Energy:                   "<<system_info.FsTotalEnergy
                         <<"\n#Total_Magnetic_Moment:          "<<system_info.FsTotalMMoment
                         <<"\n#Steps_time:                     "<<system_info.Fscycles <<std::endl;
      simple_result_vs_T <<"\n#Temperature\tVolume\tTotal_Energy\tTEnergy_error\tTotal_Moment\tTMMoment_error"<<std::endl;

      simple_result_vs_T.close();

      break;
    }

    case 2: { /*!f_print_out.opt.2 (use only to write the new data in the files. )*/

      std::ofstream simple_result_vs_T ("result_vs_T.dat", std::ios::out|std::ios::app);
      if (!simple_result_vs_T.is_open()) {
        std::wcerr << "Error abriendo \"result_vs_T.dat\" para escribir en el." << std::endl;
        abort ();
      }

      simple_result_vs_T <<system_info.FsTemperature << "\t"
                         <<system_info.FsVolume<<"\t"
                         <<system_info.FsTotalEnergy<<"\t"
                         <<0<<"\t"
                         <<system_info.FsTotalMMoment<<"\t"
                         <<0<< std::endl;

      simple_result_vs_T.close();
      break;
    }

    default: printf("f_print_out: Choice other than 1 (to create the head and remove the old files), 2 () and 3(to write the new data in the files)\n");
              break;
   }
  return 0;
}

double framework::f_liveAverage(statistic_s & system_statistic,
	                               const int step,
                                 const int opt)
{
  switch(opt){
    case 1: {/*!opt.1 option obsolete remove or change  */
      system_statistic.sSFinalAverageEnergy = 0;
      system_statistic.sSFinalAverageMMoment = 0;
      //std::cout<<"system_statistic.sSFinalAverageMMoment:"<<system_statistic.sSFinalAverageMMoment <<std::endl;
      for (int i = step; i > (step - 500); --i){
        system_statistic.sSFinalAverageEnergy += system_statistic.sSEnergy_vector[i];
        system_statistic.sSFinalAverageMMoment += system_statistic.sSMMoment_vector[i];
        //std::cout<<"i:"<<i<<"\tsystem_statistic.sSMMoment_vector[i]:"<<system_statistic.sSMMoment_vector[i]<<std::endl;
      }

      system_statistic.sSFinalAverageEnergy /= 500;
      system_statistic.sSFinalAverageMMoment /= 500;
      //std::cout<<"system_statistic.sSFinalAverageMMoment:"<<system_statistic.sSFinalAverageMMoment <<std::endl;
      return 0;
      break;
    }

    case 2: {/*!opt.2 calculate on run time the average of E/M and save each value.*/

        system_statistic.sSLiveAveMM = system_statistic.sSLiveAveMM * step + system_statistic.sSMMoment_vector[step];
        system_statistic.sSLiveAveMM /= (step + 1);
        system_statistic.sSLiveAverageMMoment[step] = system_statistic.sSLiveAveMM;

        system_statistic.sSLiveAveE = system_statistic.sSLiveAveE * step + system_statistic.sSEnergy_vector[step];
        system_statistic.sSLiveAveE /= (step + 1);
        system_statistic.sSLiveAverageEnergy[step] = system_statistic.sSLiveAveE;
        //std::cout<<system_statistic.sSEnergy_vector[step]<<std::endl;

        break;
    }
    default : std::cout<<"f_liveAverage"<<std::endl;
    break;
  }
  return 0;
}

int framework::f_detect_border(std::vector<struct_Sys_properties> & system_in,
	                             const framework::Fs_general_info & system_info)
{
  std::cout <<"\n---------- Determining the border atoms.----------"<<std::endl;
  int nn = system_info.FsNi * system_info.FsNi;
  for (int i = 0; i < system_info.FsNi; ++i){

    int Vi = int (i * nn), Vi_end = int ((i * nn) + system_info.FsNi - 1), Vf = int ((i * nn) + (nn) - system_info.FsNi), Vf_end = int ((i * nn) + (nn) - 1);

    for (int e = 0; e < (system_info.FsNi * system_info.FsNi); ++e){
      system_in[Vi + e].sBorder[0] = ((Vi + e) >= Vi && (Vi + e) <= Vi_end) ? (1) : (0 );                     //!North orientation +y
      system_in[Vi + e].sBorder[1] = ((Vi + e) >= Vf && (Vi + e) <= Vf_end) ? (1) : (0);   //!South orientation -y
      system_in[Vi + e].sBorder[2] = (!((Vi + e +1) % system_info.FsNi)) ? (1) : (0);                                 //!East orientation +x
      system_in[Vi + e].sBorder[3] = (!((Vi + e ) % system_info.FsNi)) ? (1) : (0);                                  //!West orientation -x
      system_in[Vi + e].sBorder[4] = ((Vi + e) >= ((system_info.FsNi - 1) * nn)) ? (1) : (0);  //!Up orientation +z
      system_in[Vi + e].sBorder[5] = ((Vi + e) < (nn)) ? (1) : (0);             //!Dow orientation -z
      /*------------------ to print each value  ------------------ */

        /*std::cout<<system_in[Vi + e].sBorder[0] <<"\t";
          std::cout<<system_in[Vi + e].sBorder[1] <<"\t";
          std::cout<<system_in[Vi + e].sBorder[2] <<"\t";
          std::cout<<system_in[Vi + e].sBorder[3] <<"\t";
          std::cout<<system_in[Vi + e].sBorder[4] <<"\t";
          std::cout<<system_in[Vi + e].sBorder[5] <<std::endl;
        */
    }
  }
  std::puts("Done!");
  return 0;
}

double framework::f_tow_p_distance(const std::vector<struct_Sys_properties> & system_in,
                                   const Fs_general_info & system_info,
                                   const std::vector<std::vector<int> > & neighbour,
                                   const int index,
                                   const int i)
{
	double dx = 0.0,
         dy = 0.0,
         dz = 0.0, r;
	double neighbour_position_x = system_in[neighbour[index][i]].x,
	       neighbour_position_y = system_in[neighbour[index][i]].y,
	       neighbour_position_z = system_in[neighbour[index][i]].z;

         //std::cout<<"antes del calculo\t"<<neighbour_position_x<<"\t"<<neighbour_position_y<<"\t"<<neighbour_position_z<<std::endl;

  switch (i + 1){
    case 1: {/*! to get back the correct +y distance or North distance orientation .*/
      neighbour_position_y = (system_in[index].sBorder[0]) ? (system_in[index].y + (system_info.FsL - system_in[neighbour[index][i]].y)) : ((system_in[neighbour[index][i]].y));
      break;
    }

    case 2: {/*! to get back the correct -y distance or South  */
      neighbour_position_y = (system_in[index].sBorder[1]) ? (system_info.FsL + system_in[neighbour[index][i]].y) : ((system_in[neighbour[index][i]].y));
      break;
    }

    case 3: {/*! to get back the correct +x distance or East.*/
      neighbour_position_x = (system_in[index].sBorder[2]) ? (system_in[index].x + (system_info.FsL - system_in[neighbour[index][i]].x)) : ((system_in[ neighbour[index][i]].x));
      break;
    }

    case 4: {/*! to get back the correct -x distance or West. */
      //std::cout<<"case4\n";
      neighbour_position_x = (system_in[index].sBorder[3]) ? (system_info.FsL + system_in[neighbour[index][i]].x) : ((system_in[ neighbour[index][i]].x));

      break;
    }

    case 5: {/*! to get back the correct +z distance or Up.*/
      neighbour_position_z = (system_in[index].sBorder[4]) ? (system_in[index].z + (system_info.FsL - system_in[neighbour[index][i]].z)) : ((system_in[neighbour[index][i]].z));
      //std::cout<<"neighbour_position_z:"<<neighbour_position_z<<std::endl;
      break;
    }

    case 6: {/*! to get back the correct -z distance or Down. */
      neighbour_position_z = (system_in[index].sBorder[5]) ? (system_info.FsL + system_in[neighbour[index][i]].z) : ((system_in[neighbour[index][i]].z));
      break;
    }

    default: printf("Choice other than 1, 2, 3, 4, 5 or 6. \n");
    break;
  }

  dx = float((system_in[index].x - neighbour_position_x ) * (system_in[index].x - neighbour_position_x));
  dy = float((system_in[index].y - neighbour_position_y ) * (system_in[index].y - neighbour_position_y));
  dz = float((system_in[index].z - neighbour_position_z) * (system_in[index].z - neighbour_position_z));
  r = sqrt(dx + dy + dz);

  return (float)r;
}


 double framework::f_tow_p_distance(const std::vector<struct_Sys_properties> & system_in,
                                    const struct_Sys_properties & system_tmp,
                                    const Fs_general_info & system_info,
                                    const std::vector<std::vector<int> > & neighbour,
                                    const int index,
                                    const int i)
{
  double dx = 0.0,
         dy = 0.0,
         dz = 0.0, r;
  double neighbour_position_x = system_in[neighbour[index][i]].x,
         neighbour_position_y = system_in[neighbour[index][i]].y,
         neighbour_position_z = system_in[neighbour[index][i]].z;

         //std::cout<<"antes del calculo\t"<<neighbour_position_x<<"\t"<<neighbour_position_y<<"\t"<<neighbour_position_z<<std::endl;

  switch (i + 1){
    case 1: {/*! to get back the correct +y distance or North .*/
      neighbour_position_y = (system_in[index].sBorder[0]) ? (system_tmp.y + (system_info.FsL - system_in[neighbour[index][i]].y)) : ((system_in[neighbour[index][i]].y));
      //std::cout<<"case1\n";
      break;
    }

    case 2: {/*! to get back the correct -y distance or South  */
      neighbour_position_y = (system_in[index].sBorder[1]) ? (system_info.FsL + system_in[neighbour[index][i]].y) : ((system_in[neighbour[index][i]].y));
      //std::cout<<"case2\n";
      break;
    }

    case 3: {/*! to get back the correct +x distance or East.*/
      neighbour_position_x = (system_in[index].sBorder[2]) ? (system_tmp.x + (system_info.FsL - system_in[neighbour[index][i]].x)) : ((system_in[neighbour[index][i]].x));
      break;
    }

    case 4: {/*! to get back the correct -x distance or West. */
      neighbour_position_x = (system_in[index].sBorder[3]) ? (system_info.FsL + system_in[neighbour[index][i]].x) : ((system_in[ neighbour[index][i]].x));

      break;
    }

    case 5: {/*! to get back the correct +z distance or Up.*/
      neighbour_position_z = (system_in[index].sBorder[4]) ? (system_tmp.z + (system_info.FsL - system_in[neighbour[index][i]].z)) : ((system_in[neighbour[index][i]].z));
      break;
    }

    case 6: {/*! to get back the correct -z distance or Down. */
      neighbour_position_z = (system_in[index].sBorder[5]) ? (system_info.FsL + system_in[neighbour[index][i]].z) : ((system_in[neighbour[index][i]].z));
      break;
    }

    default: printf("Choice other than 1, 2, 3, 4, 5 or 6. \n");
    break;
  }

  dx = float((system_tmp.x - neighbour_position_x ) * (system_tmp.x - neighbour_position_x));
  dy = float((system_tmp.y - neighbour_position_y ) * (system_tmp.y - neighbour_position_y));
  dz = float((system_tmp.z - neighbour_position_z) * (system_tmp.z - neighbour_position_z));

  r = sqrt(dx + dy + dz);
  //std::cout<<"------------------------------------------------"<<std::endl;
  return (float)r;
}


double framework::f_STD_desviation(Fs_general_info& system_info,
                                    statistic_s & system_statistic,
                                    const int step,
                                    const int opt)
{
  switch(opt){

    case 1:{/*f_STD_desviation:case 1. Determine the STD Error of spin and energy taking the last 500 steps.*/
      double tmp_average_spin = 0,
             tmp_average_energy = 0;

      for (int i = step; i > (step - 500); --i){
        tmp_average_spin += system_statistic.sSMMoment_vector[i];
        //std::cout<<i<<"\tsSMMoment_vector[1]"<<system_statistic.sSMMoment_vector[i]<<std::endl;
        tmp_average_energy += system_statistic.sSEnergy_vector[i];
      }

      std::cout<<"tmp_average_spin:"<<tmp_average_spin<<std::endl;

      tmp_average_spin /= 500;
      tmp_average_energy /= 500;


      for (int i = 0; i < 500; ++i){
        system_statistic.sSMoment_STD_error += (system_statistic.sSMMoment_vector[i] - tmp_average_spin) *
                                               (system_statistic.sSMMoment_vector[i] - tmp_average_spin);

        system_statistic.sSEnergy_STD_error += (system_statistic.sSEnergy_vector[i] - tmp_average_energy) *
                                               (system_statistic.sSEnergy_vector[i] - tmp_average_energy);
      }

      system_statistic.sSMoment_STD_error = sqrt( system_statistic.sSMoment_STD_error/500);
      std::cout<<"sSMoment_STD_error:"<<system_statistic.sSMoment_STD_error<<std::endl;
      system_statistic.sSEnergy_STD_error = sqrt( system_statistic.sSEnergy_STD_error/500);
      std::cout<<"sSEnergy_STD_error:"<<system_statistic.sSEnergy_STD_error<<std::endl;

      break;
    }

    case 2:{/*f_STD_desviation:case 2. Determine the live STD Error of spin and energy taking.*/

      double tmp_STD_Moment = 0, tmp_STD_Energy = 0;
      tmp_STD_Moment = (system_statistic.sSMMoment_vector[step] - system_statistic.sSLiveAverageMMoment[step]) *
                       (system_statistic.sSMMoment_vector[step] - system_statistic.sSLiveAverageMMoment[step]);

      tmp_STD_Energy = (system_statistic.sSEnergy_vector[step] - system_statistic.sSLiveAverageEnergy[step]) *
                       (system_statistic.sSEnergy_vector[step] - system_statistic.sSLiveAverageEnergy[step]);

      system_statistic.sSMoment_STD_error_vector[step] = sqrt(system_statistic.sSMoment_STD_error/(step + 1));
      system_statistic.sSEnergy_STD_error_vector[step] = sqrt(system_statistic.sSEnergy_STD_error/(step + 1));

      break;
    }
    default: std::cout<<"f_STD_desviation"<<std::endl;
  }

  return 0.0;
}

double framework::f_test_function(std::vector<struct_Sys_properties> &system_in,           /*!const.vector type structure*/
                                  Fs_general_info & system_info,
                                  s_options & system_options,                                   /*!*/
                                  const std::vector<std::vector<int> > & neighbour,               /*!2D int vector */
                                  statistic_s & system_statistic,
                                  const int index,                              /*!framework::struct_Sys_properties temporal recipient*/
                                  const int opt)
{
  switch (opt){
    case 1: {/*! f_test_function opt 1. print on the screen each atom and their neighbour.*/
      std::cout <<"----------------- neighbour --------------------"<<std::endl;
      for (int i = 0; i < system_info.FsTotalAtoms; ++i){
        std::cout<<i<<"\t";
        for (int e = 0; e < 6; ++e) std::cout<<neighbour[i][e]<<"\t";
        std::cout<<std::endl;
      }

      break;
    }

    case 2: {/*! f_test_function opt 2. print on the screen the each atom Spin.*/

      std::cout <<"----------------- lattice --------------------"<<std::endl;
      std::cout<<"Energy of spin:"<< (f_totalEnergy_totalMoment(system_in, neighbour, system_info, system_statistic, index, 5)) <<std::endl;
      for (int i = 0; i < system_info.FsTotalAtoms/system_info.FsNi; ++i){
        for (int e = 0; e < system_info.FsNi; ++e) std::cout<< system_in[e].tSpin <<"\t" ;
        std::cout<<std::endl;
      }

      break;
    }

    case 3: {/*! f_test_function opt 3. print on the screen each atom and if is on the border or not, using the code yes (1) or not (0).*/
      for (int i = 0; i < system_info.FsTotalAtoms; ++i){
        for (int e = 0; e < 6; ++e) std::cout<<i<<"\t"<<system_in[i].sBorder[e] <<"\t" ;
        std::cout<< std::endl;
      }
      break;
    }

    case 4: {/*! f_test_function opt 4. this case run only one metropolis cycle.*/
      f_metropolisMethod(system_in, neighbour, system_info, system_statistic, system_options, system_options.metropolisMethod);
      break;
    }

    case 5: {/*! f_test_function opt 5. this case run only one metropolis cycle.*/

      f_routine(system_in, neighbour, system_info, system_statistic, system_options, system_options.metropolisMethod, system_info.Fscycles);
      break;
    }

    default: printf("f_test_function.  Choice other than 1, 2 and 3\n");
    break;
  }
  return 0;
}


void framework::f_about_system_inf(const Fs_general_info & system_info,
                                   const int opt)
{
  switch(opt){
    case 1:{/*! opt 1. Init info.*/
      std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
      std::time_t now_c = std::chrono::system_clock::to_time_t(now- std::chrono::hours(24));
      std::cout <<"------------------" <<  std::put_time(std::localtime(&now_c), "%F %T")  <<" ------------------"<<std::endl;
      std::cout <<"----------- Registered the information system -----------\n"
                <<"Temperature:          " << system_info.FsTemperature<<"\n"
                <<"Pressure:             " << system_info.FsPressure<<"\n"
                <<"Volume:               " << system_info.FsVolume <<"\n"
                <<"L:                    " << system_info.FsL <<"\n"
                <<"Total_Energy:         " << system_info.FsTotalEnergy <<"\n"
                <<"Total_MMoment:        " << system_info.FsTotalMMoment <<"\n"
                <<"K_nn:                 " << system_info.Fsk_nn <<"\n"
                <<"T_min:                " << system_info.FsMinT <<"\n"
                <<"T_max:                " << system_info.FsMaxT <<"\n"
                <<"T_steps:              " << system_info.FsTsteps <<"\n"
                <<"D:                    " << system_info.FsD <<"\n"
                <<"j_isin:               " << system_info.Fsj_isin <<"\n"
                <<"Energy_Error:         " << system_info.FsTEnergy_error <<"\n"
                <<"Total_MMoment_Error:  " << system_info.FsTMMoment_error <<"\n"
                <<"R_hs:                 " << system_info.FsR_hs <<"\n"
                <<"R_ls:                 " << system_info.FsR_ls <<"\n"
                <<"Total_Atoms:          " << system_info.FsTotalAtoms <<"\n"
                <<"N_xyz:                " << system_info.FsNi <<"\n"
                <<"Cycles:               " << system_info.Fscycles <<"\n"
                <<"---------------------------------------------------------\n"
                //<< "---NOTE: the labels with (*) has been defined by the user ---"
                << std::endl;
      break;
    }

    case 2:{/*! opt 2. Final info.*/
      std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
      std::time_t now_c = std::chrono::system_clock::to_time_t(now- std::chrono::hours(24));
      std::cout <<"------------ Work has been done! -------------"<<std::endl;
      std::cout <<"------------------ " <<  std::put_time(std::localtime(&now_c), "%F %T")  <<" ------------------"<<std::endl;
      std::cout <<"----------- Registered the information system -----------\n"
                <<"Temperature:          " << system_info.FsTemperature<<"\n"
                <<"Pressure:             " << system_info.FsPressure<<"\n"
                <<"Volume:               " << system_info.FsVolume <<"\n"
                <<"L:                    " << system_info.FsL <<"\n"
                <<"Total_Energy:         " << system_info.FsTotalEnergy <<"\n"
                <<"Total_MMoment:        " << system_info.FsTotalMMoment <<"\n"
                <<"K_nn:                 " << system_info.Fsk_nn <<"\n"
                <<"T_min:                " << system_info.FsMinT <<"\n"
                <<"T_max:                " << system_info.FsMaxT <<"\n"
                <<"T_steps:              " << system_info.FsTsteps <<"\n"
                <<"D:                    " << system_info.FsD <<"\n"
                <<"j_isin:               " << system_info.Fsj_isin <<"\n"
                <<"Energy_Error:         " << system_info.FsTEnergy_error <<"\n"
                <<"Total_MMoment_Error:  " << system_info.FsTMMoment_error <<"\n"
                <<"R_hs:                 " << system_info.FsR_hs <<"\n"
                <<"R_ls:                 " << system_info.FsR_ls <<"\n"
                <<"Total_Atoms:          " << system_info.FsTotalAtoms <<"\n"
                <<"N_xyz:                " << system_info.FsNi <<"\n"
                <<"Cycles:               " << system_info.Fscycles <<"\n"
                <<"---------------------------------------------------------\n"
                //<< "---NOTE: the labels with (*) has been defined by the user ---"
                << std::endl;
      break;
    }
    default: std::cout<<""<<std::endl;
  }
}

int framework::f_creating_headfileout(const int opt)
{

  switch(opt){
    case 1:{/*opt 1. remove the VMD old file.*/

      break;
    }
    case 2:{/*opt 2. creating the head of data per temperature file and remove the old file.*/

      //std::cout <<"\n------------------- remove/create files of data  -------------------"<<std::endl;

      std::ofstream energy_by_stepFlux ("energy_by_step.dat", std::ios::out|std::ios::app);
      if (!energy_by_stepFlux.is_open()) {
        std::wcerr << "Error abriendo \"energy_by_step.dat\" para escribir en el." << std::endl;
        abort ();
      }
      energy_by_stepFlux <<"#T\tstep \tL_Energy\tEnergy_LAverage\tL_T_M_Moment\tT_M_Moment_LAverage\tEnergy\tEnergy_STDE\tMagnetic_Moment\tMagnetic_Moment_STDE"<<std::endl;
      energy_by_stepFlux.close();

      break;
    }
    case 3:{/*opt 3. creating of a big data file*/

      break;
    }

    default: std::cout<<"Try with another option."<<std::endl;
  }

  return 0;
}


void framework::f_remove_olf_files()
{

  if(remove("isin.dat") != 0) std::perror("Error deleting file: isin.dat");
    else std::puts("The ´isin.dat´ file was successfully deleted.");

  if(remove("vmd_file.xyz") != 0) std::perror("Error deleting file: vmd_file.xyz");
    else std::puts("The ´vmd_file.xyz´ file was successfully deleted.");

  if(remove("system.dat") != 0) std::perror("Error deleting file: system.dat");
    else std::puts("Old file successfully deleted");

  if(remove("system_info.dat") != 0) std::perror("Error deleting file: system_info.dat");
    else std::puts("Old file successfully deleted");

  if(remove("system_info_plot.dat") != 0) std::perror("Error deleting file: system_info_plot.dat");
    else std::puts("Old file successfully deleted");

  if(remove("result_vs_T.dat") != 0) std::perror("Error deleting file: result_vs_T.dat");
    else std::puts("Old file successfully deleted");

  if(remove("time_xyz.dat") != 0) std::perror("Error deleting file: time_xyz.dat");
    else std::puts("Old file successfully deleted");

  if(remove("Data.dat") != 0) std::perror("Error deleting file: time_xyz.dat");
    else std::puts("Old file successfully deleted");

  if(remove("liveData.dat") != 0) std::perror("Error deleting file: liveData.dat");
    else std::puts("Old file successfully deleted");

  if(remove("energy_by_step.dat") != 0) std::perror("Error deleting file: energy_by_step.dat");
    else std::puts("Old file successfully deleted");
}



/*!******************************************************************************************/
/*! f_temperature_file_print(): public function to print in the file out called "Data.dat"  */
/*!                   the final result after isotherm stabilization                         */
/*!******************************************************************************************/

int framework::f_temperature_file_print(const std::vector<struct_Sys_properties> & system_in,
                                        const Fs_general_info & system_info,
                                        const statistic_s & system_statistic,
                                        const int opt)
{
  switch (opt){
    case 1: {/*! f_temperature_file_print opt 1. print Head.(bigdata)*/

      //std::cout <<"\n------------------- remove/create files of data  -------------------"<<std::endl;

      std::ofstream FData ("Data.dat", std::ios::out|std::ios::app);
      if (!FData.is_open()) {
        std::wcerr << "Error abriendo \"Data.dat\" para escribir en el." << std::endl;
        abort ();
      }

      FData <<"#T\tEnergy_Average\tEnergy_Error\tM_Moment_Average\tM_Moment_Error"<<std::endl;

      FData.close();

      break;
    }

    case 2: {/*! f_temperature_file_print opt 2. print the data. (bigdata) */

      std::ofstream dataFlux ("Data.dat", std::ios::out|std::ios::app);
      if (!dataFlux.is_open()) {
        std::wcerr << "Error abriendo \"Data.dat\" para escribir en el." << std::endl;
        abort ();
      }

      dataFlux <<system_info.FsTemperature<<"\t"
               <<system_statistic.sSFinalAverageEnergy <<"\t"
               <<system_statistic.sSEnergy_STD_error <<"\t"
               <<system_statistic.sSFinalAverageMMoment <<"\t"
               <<system_statistic.sSMoment_STD_error <<std::endl;
      //dataFlux <<std::endl;

      dataFlux.close();


      break;
    }

    default: printf("f_temperature_file_print.  Choice other than 1, 2 \n");
    break;
  }
  return 0;
}


/*!*********************************************************************************/
/*! f_xyz_print(): public function to print in the file out called "time_xyz.dat"  */
/*!                   the result at the end of one MC method.                      */
/*!*********************************************************************************/

int framework::f_xyz_print(const std::vector<struct_Sys_properties> & system_in,
                           const Fs_general_info & system_info,
                           const statistic_s & system_statistic,
                           const int step,
                           const int opt)
{
  switch (opt){
    case 1: {/*! f_xyz_print opt 1. print only Head.*/

      std::cout <<"\n------------------- create files of data  -------------------"<<std::endl;
      std::ofstream xyzFlux ("time_xyz.dat", std::ios::out|std::ios::app);
      if (!xyzFlux.is_open()) {
        std::wcerr << "Error abriendo \"time_xyz.dat\" para escribir en el." << std::endl;
        abort ();
      }

      xyzFlux <<"#step\tn\tT\tx\ty\tz\tR\tS"<<std::endl;
      xyzFlux.close();

      break;
    }

    case 2: {/*! f_xyz_print opt 2. print the data without head */

      std::ofstream xyzFlux ("time_xyz.dat", std::ios::out|std::ios::app);
      if (!xyzFlux.is_open()) {
        std::wcerr << "Error abriendo \"time_xyz.dat\" para escribir en el." << std::endl;
        abort ();
      }

      for (int i = 0; i < system_info.FsTotalAtoms; ++i)  xyzFlux <<step<<"\t"
                                                                  <<i<<"\t"
                                                                  << system_info.FsTemperature << "\t"
                                                                  << system_in[i].x <<"\t"
                                                                  << system_in[i].y <<"\t"
                                                                  << system_in[i].z <<"\t"
                                                                  << system_in[i].R <<"\t"
                                                                  << system_in[i].tSpin << std::endl;
      xyzFlux <<std::endl;
      xyzFlux.close();
      break;
    }

    default: printf("f_xyz_print.  Choice other than 1, 2 and 3\n");
    break;
  }

  return 0;
}


/*!****************************************************************************************/
/*! f_cycle_file_print(): public function to print in the file out called "liveData.dat"  */
/*!                   the result of each iteration cycle of MC.                           */
/*!****************************************************************************************/

int framework::f_cycle_file_print(const std::vector<struct_Sys_properties> & system_in,
                                  const Fs_general_info & system_info,
                                  const statistic_s & system_statistic,
                                  const s_options & system_options,
                                  const int opt,
                                  const int step)
{
  switch (opt){
    case 1: {/*! f_cycle_file_print opt 1. print Head.*/

      //std::cout <<"\n------------------- remove/create files of data  -------------------"<<std::endl;

      std::ofstream liveData("liveData.dat", std::ios::out|std::ios::app);
      if (!liveData.is_open()) {
        std::wcerr << "Error abriendo \"liveData.dat\" para escribir en el." << std::endl;
        abort ();
      }

      std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
      std::time_t now_c = std::chrono::system_clock::to_time_t(now- std::chrono::hours(24));
      liveData  <<"#------------------" <<  std::put_time(std::localtime(&now_c), "%F %T")  <<" ------------------"<<std::endl;
      liveData  <<"#----------- Registered the information system -----------\n"
                <<"#Temperature:          " << system_info.FsTemperature<<"\n"
                <<"#Pressure:             " << system_info.FsPressure<<"\n"
                <<"#Volume:               " << system_info.FsVolume <<"\n"
                <<"#L:                    " << system_info.FsL <<"\n"
                <<"#Total_Energy:         " << system_info.FsTotalEnergy <<"\n"
                <<"#Total_MMoment:        " << system_info.FsTotalMMoment <<"\n"
                <<"#K_nn:                 " << system_info.Fsk_nn <<"\n"
                <<"#T_min:                " << system_info.FsMinT <<"\n"
                <<"#T_max:                " << system_info.FsMaxT <<"\n"
                <<"#T_steps:              " << system_info.FsTsteps <<"\n"
                <<"#D:                    " << system_info.FsD <<"\n"
                <<"#j_isin:               " << system_info.Fsj_isin <<"\n"
                <<"#Energy_Error:         " << system_info.FsTEnergy_error <<"\n"
                <<"#Total_MMoment_Error:  " << system_info.FsTMMoment_error <<"\n"
                <<"#R_hs:                 " << system_info.FsR_hs <<"\n"
                <<"#R_ls:                 " << system_info.FsR_ls <<"\n"
                <<"#Total_Atoms:          " << system_info.FsTotalAtoms <<"\n"
                <<"#N_xyz:                " << system_info.FsNi <<"\n"
                <<"#Cycles:               " << system_info.Fscycles <<"\n"
                <<"#Routine_case:         " << system_options.routine<<"\n"
                <<"#Metropolis_case:      " << system_options.metropolisMethod<<"\n"
                <<"#Switch_Energy_case:   " << system_options.swithEnergy<<"\n"
                <<"#---------------------------------------------------------\n"
                << std::endl;


      liveData <<"#step\tTemperature\tL_Energy\tEnergy_LAverage\tEnergy_Error\tL_T_M_Moment\tT_M_Moment_LAverage\tT_M_Moment_Error"<<std::endl;

      liveData.close();

      break;
    }

    case 2: {/*! f_cycle_file_print opt 2. print the data.*/

      std::ofstream liveData ("liveData.dat", std::ios::out|std::ios::app);
      if (!liveData.is_open()) {
        std::wcerr << "Error abriendo \"liveData.dat\" para escribir en el." << std::endl;
        abort ();
      }

      liveData<<step<<"\t"
              <<system_info.FsTemperature<<"\t"
              <<system_statistic.sSEnergy_vector[step] <<"\t"
              <<system_statistic.sSLiveAverageEnergy[step] <<"\t"
              <<system_statistic.sSEnergy_STD_error_vector[step]<<"\t"
              <<system_statistic.sSMMoment_vector[step]<<"\t"
              <<system_statistic.sSLiveAverageMMoment[step]<<"\t"
              <<system_statistic.sSMoment_STD_error_vector[step]<<std::endl;
      liveData.close();


      break;
    }

    case 3: {/*! f_cycle_file_print opt 3. to print the last line */

      std::ofstream liveData ("liveData.dat", std::ios::out|std::ios::app);
      if (!liveData.is_open()) {
        std::wcerr << "Error abriendo \"liveData.dat\" para escribir en el." << std::endl;
        abort ();
      }

      liveData <<"LC"<<"\t"
               <<system_info.FsTemperature<<"\t"
               <<system_statistic.sSEnergy_vector[step] <<"\t"
               <<system_statistic.sSLiveAverageEnergy[step] <<"\t"
               <<system_statistic.sSEnergy_STD_error_vector[step]<<"\t"
               <<system_statistic.sSMMoment_vector[step]<<"\t"
               <<system_statistic.sSLiveAverageMMoment[step]<<"\t"
               <<system_statistic.sSMoment_STD_error_vector[step]<<"\t"
               <<system_statistic.sSFinalAverageEnergy<<"\t"
               <<system_statistic.sSEnergy_STD_error<<"\t"
               <<system_statistic.sSFinalAverageMMoment<<"\t"
               <<system_statistic.sSMoment_STD_error<<"\t"<<std::endl;
      std::cout<<"Moment:"<<system_statistic.sSFinalAverageMMoment<<std::endl;



      liveData <<std::endl;
      liveData.close();
      break;
    }

    default: printf("f_cycle_file_print: Choice other than 1, 2 or 3. %d \n", opt);
    break;
  }
  return 0;
}


/*******************************************/
/* f_print_fileout: public function to     */
/*      print the neighbour vector         */
/*******************************************/

int framework::f_print_fileout (const std::vector<struct_Sys_properties> &system_in,
                              const std::vector<std::vector<int> > &neighbour,
                              Fs_general_info & system_info,
                              const statistic_s & system_statistic,
                              const int opt,
                              const int cycle)
{
  switch (opt){

    case 1: { /*!opt.1 to print on the screen the neighbour of [i] atom.*/
      std::cout<<"------------- printing on the screen the neighbour of [i] atom. -------------\n";
      for (int i = 0; i < neighbour.size(); ++i){
        std::cout <<"i"<< i << "\t"
                  <<"N"<< neighbour[i][0]<< "\t"
                  <<"S"<< neighbour[i][1]<< "\t"
                  <<"E"<< neighbour[i][2]<< "\t"
                  <<"W"<< neighbour[i][3]<< "\t"
                  <<"U"<< neighbour[i][4]<< "\t"
                  <<"D"<< neighbour[i][5]<< "\t"<< std::endl;
      }
      break;
    }

    case 2: { /*!opt.2 this case shall remove/create three files of data.*/

      std::cout <<"\n------------------- remove/create files of data  -------------------"<<std::endl;

      std::ofstream outflux ("system.dat", std::ios::out|std::ios::app);
      if (!outflux.is_open()) {
        std::wcerr << "Error abriendo \"system.dat\" para escribir en el." << std::endl;
        abort ();
      }

      outflux <<"#step\tn\tx\ty\tz\tS\tR"<<std::endl;

      outflux.close();

      /********************************* to print system_info ****************************************/



      std::ofstream p_systeminfo ("system_info.dat", std::ios::out|std::ios::app);
      if (!p_systeminfo.is_open()) {
        std::wcerr << "Error abriendo \"system_info.dat\" para escribir en el." << std::endl;
        abort ();
      }
      p_systeminfo <<"Cycle:             "<< cycle<<"\n"
                   <<"Total_Atoms:       "<< system_info.FsTotalAtoms<<"\n"
                   <<"n per side:        "<< system_info.FsNi<<"\n"
                   <<"Temperature:       "<< system_info.FsTemperature<<"\n"
                   <<"Pressure:          "<< system_info.FsPressure<<"\n"
                   <<"Volume:            "<< system_info.FsVolume<<"\n"
                   <<"L:                 "<< system_info.FsL<<"\n"
                   <<"Total_Energy:      "<< system_info.FsTotalEnergy<<"\n"
                   <<"Tota_MMoment:      "<< system_info.FsTotalMMoment <<std::endl;

      p_systeminfo.close();

      /********************************* to print system_info_plot ****************************************/



      std::ofstream outflux_syst_inf_plot ("system_info_plot.dat", std::ios::out|std::ios::app);
      if (!outflux_syst_inf_plot.is_open()) {
        std::wcerr << "Error abriendo \"system_info_plot.dat\" para escribir en el." << std::endl;
        abort ();
      }
      outflux_syst_inf_plot <<"#n_cyle\tTotal_Atoms\tsqurt(N)\tTemperature\tPressure\tL\tVolume\tTotal_Energy\tAverage_Energy\tTotal_Spin\tAverage_MMoment"<<std::endl;

      outflux_syst_inf_plot.close();

      //std::puts("Done!");
      break;
    }

    case 3: { /*!opt.3 this case print in file the data.*/

      std::ofstream outflux ("system.dat", std::ios::out|std::ios::app);
      if (!outflux.is_open()) {
        std::wcerr << "Error abriendo \"outputfile.log\" para escribir en el." << std::endl;
        abort ();
      }

      for (int i = 0; i < (system_info.FsTotalAtoms); ++i)  outflux << cycle << "\t"
                                                                    << i << "\t"
                                                                    << system_in[i].x << "\t"
                                                                    << system_in[i].y << "\t"
                                                                    << system_in[i].z <<"\t"
                                                                    << system_in[i].tSpin <<"\t"
                                                                    << system_in[i].R<< std::endl;
      outflux.close();

      /********************************* to print system_info ****************************************/

      std::ofstream outflux_systeminfo ("system_info.dat", std::ios::out|std::ios::app);
      if (!outflux_systeminfo.is_open()) {
        std::wcerr << "Error abriendo \"system_info.dat\" para escribir en el." << std::endl;
        abort ();
      }
      outflux_systeminfo <<std::endl;
      outflux_systeminfo <<"Cycle:             "<< cycle<<"\n"
                         <<"Total_Atoms:       "<< system_info.FsTotalAtoms<<"\n"
                         <<"n per side:        "<< system_info.FsNi<<"\n"
                         <<"Temperature:       "<< system_info.FsTemperature<<"\n"
                         <<"Pressure:          "<< system_info.FsPressure<<"\n"
                         <<"Volume:            "<< system_info.FsVolume<<"\n"
                         <<"L:                 "<< system_info.FsL<<"\n"
                         <<"Total_Energy:      "<< system_info.FsTotalEnergy<<"\n"
                         <<"Tota_MMoment:      "<< system_info.FsTotalMMoment <<std::endl;
      outflux_systeminfo.close();

      /********************************* to print system_info_plot ****************************************/
      std::ofstream outflux_syst_inf_plot ("system_info_plot.dat", std::ios::out|std::ios::app);
      if (!outflux_syst_inf_plot.is_open()) {
        std::wcerr << "Error abriendo \"system_info_plot.dat\" para escribir en el." << std::endl;
        abort ();
      }
      //outflux_syst_inf_plot <<"#ncyle \t Total_Atoms \t n per side \t Temperature \t Pressure\t Volume\t L\t Total_Energy \t"<<std::endl;
      outflux_syst_inf_plot << cycle<<"\t"
                            << system_info.FsTotalAtoms<<"\t"
                            << system_info.FsNi<<"\t"
                            << system_info.FsTemperature<<"\t"
                            << system_info.FsPressure<<"\t"
                            << system_info.FsL<<"\t"
                            << system_info.FsVolume<<"\t"
                            << system_info.FsTotalEnergy<<"\t"
                            << system_statistic.sSLiveAverageEnergy[cycle]<<"\t"
                            << system_info.FsTotalMMoment
                            << system_statistic.sSLiveAverageMMoment[cycle]<<"\t"
                            <<std::endl;
      outflux_syst_inf_plot.close();

      //std::puts("Done!");
      break;
    }

    default: printf("Choice other than 1, 2 and 3\n");
    break;
  }
  return 0;
}

int framework::f_write_vector(std::vector<framework::struct_Sys_properties> &write_vector,
                            const int n)
{
  /*write_vector.resize(n * n * n);

  for (int i = 0; i < (n * n * n); ++i){
    write_vector[i].x = i;
    write_vector[i].y = i;
    write_vector[i].z = i;
    write_vector[i].tSpin=1;
  } */
  return 0;
}

int framework::f_read_vector(const std::vector<struct_Sys_properties> &read_vector,
                             const int N)
{
  std::ofstream outflux ("write_vector.dat", std::ios::out|std::ios::app);
    if (!outflux.is_open()) {
      std::wcerr << "Error abriendo \"write_vector.dat\" para escribir en el." << std::endl;
      abort ();
  }

  for (int i = 0; i < N; ++i){
    outflux << i << "\t"<< read_vector[i].x << "\t" << read_vector[i].y << "\t" << read_vector[i].z <<"\t"<< read_vector[i].R<<std::endl;
  }
  outflux <<"\n"<<std::endl;
  return 0;
}

/************************************/
/*                                  */
/*   AREA OF AEUXIALIAR FUNCTIONS   */
/* MOVE AS POSIBLE TO EXTERNAL FILE */
/*                                  */
/************************************/

int neighbourDefine (const std::vector<framework::struct_Sys_properties> & system_in,
	                   std::vector<std::vector<int> > & neighbour,
	                   const int n)
{
  std::cout <<"\n---------- Determining the neighbour of each atoms.----------"<<std::endl;
  for (int i = 0; i < n; ++i){
    int Vi = int (i * n * n), Vi_end = int ((i * n * n) + n - 1), Vf = int ((i * n * n) + (n * n) - n), Vf_end = int ((i * n * n) + (n * n) - 1);

    /*! this loop find the nearest neighbour and write them in the neighbour_matrix*/
  	for (int e = 0; e < (n * n); ++e){
      neighbour[Vi + e][0] = ((Vi + e) >= Vi && (Vi + e) <= Vi_end) ? (Vf + e) : ((Vi + e) - n);                     //!North orientation +y
      neighbour[Vi + e][1] = ((Vi + e) >= Vf && (Vi + e) <= Vf_end) ? ((Vi + e) - (n * (n - 1))) : ((Vi + e) + n);   //!South orientation -y
      neighbour[Vi + e][2] = (!((Vi + e +1) % n)) ? (Vi + e - n + 1) : (Vi + e + 1);                                 //!East orientation +x
      neighbour[Vi + e][3] = (!((Vi + e ) % n)) ? (Vi + e + n -  1) : (Vi + e - 1);                                  //!West orientation -x
      neighbour[Vi + e][4] = ((Vi + e) >= ((n - 1) * n * n)) ? ((Vi + e) - ((n - 1) * n * n)) : ((Vi + e) + n * n);  //!Up orientation +z
      neighbour[Vi + e][5] = ((Vi + e) < (n * n)) ? ((Vi + e) + ((n - 1) * n * n)) : ((Vi + e) - n * n);             //!Dow orientation -z
  	}
  }
  std::puts("Done!");
  return 0;
}

int randomSamplePosition(std::vector<int> &work_vector,
                         const int P) //! esta funcion esta bien describirla luego
{
  for (int i = 0; i < P * P * P; ++i){
    int swap, swap_two;
    swap = work_vector[i];
    swap_two = int_randomNumber_generator(i, (P * P * P) - 1);
    work_vector[i] = work_vector [swap_two];
    work_vector[swap_two] = swap;
  }
  return 0;
}

int int_randomNumber_generator(const int min/*!min*/,
                               const int max/*!max*/) /*!to generate a real number between <const min:const max>*/
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis_int(min, max);
  return dis_int(gen);
}

double real_randomNumber_generator(const int min/*!min*/,
                                   const int max/*!max*/)/*! to generate a real number between <const min:const max> */
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis_real(min, max);
  return dis_real(gen);
}

/*!************ corner_f_lattice::PRINT THE CORNERS OF EACH SLICE OF THE LATTICE ****************/
int corner_f_lattice (const std::vector<framework::struct_Sys_properties> & system_in,
	                    std::vector<int> &corner,
	                    const int n)
{
  std::cout <<"\n---------- Determining the corners number of each slice of the system. ----------"<<std::endl;
  corner[0]=(0 * n * n);                            //ground. 000
  corner[1]=int ((0 * n * n) + n - 1);              //ground. 010
  corner[2]=int ((0 * n * n) + (n * n) - n);        //ground. 100
  corner[3]=int ((0 * n * n) + (n * n) - 1);        //ground. 110
  corner[4]=((n -1) * n * n);                       //ground. 001
  corner[5]=int (((n - 1) * n * n) + n - 1);        //ground. 011
  corner[6]=int (((n - 1) * n * n) + (n * n) - n);  //ground. 101
  corner[7]=int (((n - 1) * n * n) + (n * n) - 1);  //ground. 111
  std::puts("Done!");
  return 0;
}

/*!************ VMD_EXPORT_F::FUNCTION TO EXPORT IN VMD FORMAT ****************/
int VMD_EXPORT_F(const std::vector<framework::struct_Sys_properties> & system_in,
	                const int VMDTAtoms,
	                const int opt)
{
  switch (opt){
    case 1: { /*!VMD_EXPORT_F.opt.1. Call only to remove old file.*/

      if(remove("vmd_file.xyz") != 0) std::perror("");
       else std::puts("");
      break;
    }

    case 2: { /*!VMD_EXPORT_F.opt.2. Call to write in vmd xyz format.*/
      //std::cout <<"\n---------- Exporting vmd file!. ----------"<<std::endl;
      std::ofstream VMDflux ("vmd_file.xyz", std::ios::out|std::ios::app);
      if (!VMDflux.is_open()) {
        std::wcerr << "Error abriendo \"vmd_file.dat\"." << std::endl;
        abort ();
      }
      VMDflux << VMDTAtoms <<"\n"
              << "#Coment"<<std::endl;
      for (int i = 0; i < VMDTAtoms; ++i) VMDflux << "M\t" <<system_in[i].x << "\t"<< system_in[i].y << "\t"<< system_in[i].z <<"\t"<< system_in[i].tSpin <<std::endl;

      break;
    }

    default: printf("VMD_EXPORT_F. Choice other than 1 to remove old files or 2 to create the vmd files. \n");
    break;
   }
  return 0;
}

int parse_commandline(int argc,
                      char **argv,
                      framework::Fs_general_info & system_info,
                      framework::s_options & system_options)
{

  int c;
  int itmp;
  float ulitmp;
  double dtmp;
  char filename[256];

  while (1){

    int this_option_optind = optind ? optind : 1;
    int option_index = 0;

    static struct option long_options[] =
    {
      {"dimension",  required_argument, 0, 'd'},             /*! amount of atoms in one direction to do d*d*d = N inside.*/
      {"inputfile",  required_argument, 0, 'F'},             /*! read from input file all about the new job.*/
      {"time_rate_screen",  required_argument, 0, 's'},      /*! specifies the output rate to print on screen.*/
      {"time_rate_file",  required_argument, 0, 'a'},        /*! specifies the output rate to print in output file.*/
      {"metropolis_option",  required_argument, 0, 'm'},     /*! specifies the metropolis option to run, default opt is 2.*/
      {"mc_cycles",  required_argument, 0, 'c'},             /*! amount Monte-Carlo cycles.*/
      {"rutine_option",  required_argument, 0, 'C'},         /*! select the routine case.*/
      {"energy_method",  required_argument, 0, 'e'},         /*! select the energy case routine.*/
      {"T", required_argument, 0, 'T'},                      /*! Temperature to do an isothermal run.*/
      {"xyz_optimization", 0, 0, 'G'},                      /*! Temperature to do an isothermal run.*/
      {"spin_optimization", 0, 0, 'S'},                      /*! Temperature to do an isothermal run.*/
      {"initial_temperature", required_argument, 0, 'I'},    /*! initial temperature (0.1).*/
      {"temperature_step", required_argument, 0, 'p'},       /*! temperature step (0.2)..*/
      {"final_temperature", required_argument, 0, 'M'},      /*! final temperature(10.0).*/
      {"isin_model", 0, 0, 'i'},                             /*! to run the Isin-like model.*/
      {"j_isin",    required_argument, 0, 'j'},              /*! constant value for j_isisn, representing the energy between tow atoms spin.*/
      {"k_spin",    required_argument, 0, 'k'},              /*! constant value for spring constant k.*/
      {"energy_D",    required_argument, 0, 'D'},            /*! D represent the energy gap between two spin state.*/
      {"help", 0, 0, 'h'},                                   /*! helps option.*/
      {"version", 0, 0, 'v'},                                /*! to print the version of the program.*/
      {"?", 0, 0, '?'},                                      /*! call help function!.*/
      {0, 0, 0, 0}
    };



    c = getopt_long (argc, argv, "hpiGSd:c:C:D:T:j:k:F:s:a:m:I:p:M:e:", long_options, &option_index);

    if (c == -1) /*! c will be -1 after finish the last argument, and this sentence will take out the while.*/
      break;

    switch (c) {

      case 'd':{ /*!case d: amount of atoms in one direction to do d*d*d = N inside.*/
        if(sscanf(optarg, "%d", &itmp) != 1){
          fprintf(stderr,"Error parsing -d option\n");
          help();
          exit(-1);
        }
        else {
          system_info.FsNi = itmp;
          system_info.FsTotalAtoms = system_info.FsNi * system_info.FsNi * system_info.FsNi;
        }
        break;
      }

      case 'c':{ /*!case c: amount of atoms in one direction to do d*d*d = N inside.*/
        if(sscanf(optarg, "%d", &itmp) != 1){
          fprintf(stderr,"Error parsing -c option\n");
          help();
          exit(-1);
        }
        else system_info.Fscycles = itmp;

        break;
      }

      case 's':{ /*!case s: screen rate.*/
        if(sscanf(optarg, "%d", &itmp) != 1){
          fprintf(stderr,"Error parsing -s option\n");
          help();
          exit(-1);
        }
        else system_info.FsScreenrate = itmp;

        break;
      }

      case 'G':{ /*!case G: geometry optimization.*/
        system_options.metropolisMethod = 5;
        system_options.routine = 5;
        system_options.totalEnergy_totalMoment = 3;
        break;
      }

       case 'S':{ /*!case S: total spin optimization energy.*/


        break;
      }

      case 'C':{ /*!case C: select the routine case.*/
        if(sscanf(optarg, "%d", &itmp) != 1){
          fprintf(stderr,"Error parsing -C option\n");
          help();
          exit(-1);
        }
        else system_options.routine = itmp;

        break;
      }

      case 'm':{ /*!case s: metropolis case.*/
        if(sscanf(optarg, "%d", &itmp) != 1){
          fprintf(stderr,"Error parsing -m option\n");
          help();
          exit(-1);
        }
        else system_options.metropolisMethod = itmp;

        break;
      }

      case 'a':{ /*!case s: rate to print file out.*/
        if(sscanf(optarg, "%d", &itmp) != 1){
          fprintf(stderr,"Error parsing -a option\n");
          help();
          exit(-1);
        }
        else system_info.FsFilerate = itmp;

        break;
      }

      case 'T':{ /*! case T: Temperature to do an isothermal run.*/

        if(sscanf(optarg, "%lf", &dtmp) != 1){
          fprintf(stderr,"Error parsing -T option\n");
          help();
          exit(-1);
        }
        else system_info.FsTemperature = dtmp;
        break;
      }

      case 'I':{ /*! case I: Initial temperature Temperature.*/

        if(sscanf(optarg, "%lf", &dtmp) != 1){
          fprintf(stderr,"Error parsing -I option\n");
          help();
          exit(-1);
        }
        else system_info.FsMinT = dtmp;
        break;
      }

      case 'p':{ /*! case p: Temperature to do an isothermal run.*/

        if(sscanf(optarg, "%lf", &dtmp) != 1){
          fprintf(stderr,"Error parsing -p option\n");
          help();
          exit(-1);
        }
        else system_info.FsTsteps = dtmp;
        break;
      }

      case 'M':{ /*! case N: Temperature to do an isothermal run.*/

        if(sscanf(optarg, "%lf", &dtmp) != 1){
          fprintf(stderr,"Error parsing -N option\n");
          help();
          exit(-1);
        }
        else system_info.FsMaxT = dtmp;
        break;
      }

      case 'e':{ /*!case e: enrgy method routine case.*/
        if(sscanf(optarg, "%d", &itmp) != 1){
          fprintf(stderr,"Error parsing -e option\n");
          help();
          exit(-1);
        }
        else system_options.totalEnergy_totalMoment = itmp;

        break;
      }

      case 'j':{ /*! case j: constant value for j_isisn, representing the energy between tow atoms spin.*/
        if(sscanf(optarg, "%lf", &dtmp) != 1) {
          fprintf(stderr,"Error parsing -j option\n");
          help();
          exit(-1);
        }
        else system_info.Fsj_isin = dtmp;
        break;
      }

      case 'D':{ /*! case D: .*/
        if(sscanf(optarg, "%lf", &dtmp) != 1) {
          fprintf(stderr,"Error parsing -D option\n");
          help();
          exit(-1);
        }
        else system_info.FsD = dtmp;
        break;
      }

      case 'F':{ /*! case j: constant value for j_isisn, representing the energy between tow atoms spin.*/
        if(sscanf(optarg, "%d", &itmp) != 1) {
          fprintf(stderr,"Error parsing -K option\n"); //WARNING al escribir esta funcion tener en cuenta que los otros parametros pueden sobreescribir los parametros introducidos por esta opcion
          help();
          exit(-1);
        }
        else {
          system_options.open_parammeters_file = itmp;
        }
        break;
      }

      case 'k':{ /*! case k: constant value for spring constant k.*/
        if(sscanf(optarg, "%lf", &dtmp) != 1) {
          fprintf(stderr,"Error parsing -k option\n");
          help();
          exit(-1);
        }
        else system_info.Fsk_nn = dtmp;
        break;
      }

      case 'i':{ /*case h: to call the help function..*/
        system_options.b_temperatureCycle =2;
        system_options.metropolisMethod = 3;
        system_options.routine = 4;
        system_options.totalEnergy_totalMoment = 3;
        break;
      }

      case 'h':{ /*case h: to call the help function..*/
        help();
        exit(-1);
        break;
      }

      case '?':{ /*case h: to call the help function..*/
        help();
        exit(-1);
        break;
      }

      case 'v':{
        std::cout<<"**********************************************************************************"<<std::endl;
        std::cout<<"*   You can check for more information about the project the bitbuket wiki page. *"<<std::endl;
        std::cout<<"*        https://bitbucket.org/alejandro_rodriguez__/project_sco/wiki/Home       *"<<std::endl;
        std::cout<<"**********************************************************************************"<<std::endl;
        exit(-1);
        break;
      }

      default:/*Si el argumento no machea con ninguna opción conocida, debería ser un error en los parámetros...*/
        printf ("Try './vector --help' for more information.\n");

        exit(-1);
    }
  }
  /*Si siguen apareciendo cosas que no son argumentos, se imprimen hasta que se acaben...*/
  if (optind < argc) {
    printf ("No son opciones pero estan en ARGV: ");
    while (optind < argc) printf ("%s ", argv[optind++]);
    printf ("\n");
    exit(-1);
  }
  return 0;
}

void help()
{
    std::cout<<"Usage:\t vector [options]..."<<std::endl;
    std::cout<<"Simple project to study the spin-crossover on molecular magnets.\n"<<std::endl;
    std::cout<<"  -d, --dimension                      dimension in 1D of the system, where N = d * d * d. Note: d should be integer."<<std::endl;
    std::cout<<"  -c, --mc_cycle                       number Monte-Carlos steps. Note: c should be integer."<<std::endl;
    std::cout<<"  -F, --inputfile                      input file. Note: if the input_file option is passed, all the parameter will be avoid."<<std::endl;
    std::cout<<"  -T, --T                              Temperature."<<std::endl;
    std::cout<<"  -I, --initial_temperature            initial temperature (0.1)."<<std::endl;
    std::cout<<"  -p, --temperature_step               temperature step (0.2)."<<std::endl;
    std::cout<<"  -M, --final_temperature              final temperature(10.0)."<<std::endl;
    std::cout<<"  -i, --isin_model                     execute the program on the Isin-model 3D mode."<<std::endl;
    std::cout<<"  -G, --xyz_optimization               geometry optimization mode."<<std::endl;
    std::cout<<"  -S, --spin_optimization              spin geometry energy optimization mode."<<std::endl;
    std::cout<<"  -j, --j_isin                         constant of exchange interaction."<<std::endl;
    std::cout<<"  -m, --metropolis_option              to select the metropolis case (1, 2, 3, 4 or 5)."<<std::endl;
    std::cout<<"  -C, --rutine_option                  select the routine case between 5 cases(by default 3)."<<std::endl;
    std::cout<<"  -e, --energy_method                  select the energy method routine. Have 5 cases (by default 3)."<<std::endl;
    std::cout<<"  -k, --k_spring                       constant of spring."<<std::endl;
    std::cout<<"  -D, --energy_D                       D represent the energy gap between two spin state."<<std::endl;
    std::cout<<"  -s, --time_rate_screen               specifies the output rate to print on screen.(500 by default)"<<std::endl;
    std::cout<<"  -a, --time_rate_file                 specifies the output rate to print in output file.(10 by default)"<<std::endl;
    std::cout<<"  -h, --help,                          display this help and exit (or -?)."<<std::endl;
    std::cout<<"  -v, --version                        display the version and exit.\n\n\n"<<std::endl;

    std::cout<<"  Example:"<<std::endl;
    std::cout<<"  vector -d 15 -c 10000                      defining only the size and cycles number. In this case the temperature will be defined by the default at 1K."<<std::endl;
    std::cout<<"  vector -d 15 -c 1000 -T 2.0 -j 20          defining the size,cycles number, temperature and exchange interaction constant.\n"<<std::endl;

    std::cout<<"******************************************************************************"<<std::endl;
    std::cout<<"*   You can consul the more information about the project in bitbuket wiki   *"<<std::endl;
    std::cout<<"*     https://bitbucket.org/alejandro_rodriguez__/project_sco/wiki/Home      *"<<std::endl;
    std::cout<<"******************************************************************************"<<std::endl;

    return;
}
