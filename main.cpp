#include <framework.h>
#include <m_lattice.h>
#include <map>
#include <iterator>

int read_inputfile(std::string,
                   framework::Fs_general_info &,
                   framework::s_options &);


int main(int argc, char **argv)
{
  int n = 5;


  framework ob;

  std::vector<framework::struct_Sys_properties> mesh;
  std::vector<std::vector<int> > mesh_neighbour;
  framework::Fs_general_info mesh_system_info;
  framework::statistic_s mesh_statistic_info;
  framework::s_options mesh_options;

  //parse_commandline(argc, argv, mesh_system_info, mesh_options);

  ob.f_define_system_info(mesh_system_info, mesh_statistic_info, mesh_options);
  //ob.f_define_system_info(mesh_system_info, mesh_statistic_info, mesh_options, 1.0/*! temperature*/, 2500/*!cyclesN*/, n);

  parse_commandline(argc, argv, mesh_system_info, mesh_options);

  ob.create_lattice(mesh, mesh_neighbour, mesh_system_info, mesh_system_info.FsR_ls);
  std::cout << "mesh_system_info.Fsdelta:" <<mesh_system_info.Fsdelta <<std::endl;
  neighbourDefine(mesh, mesh_neighbour, mesh_system_info.FsNi);
  ob.f_detect_border(mesh, mesh_system_info);
  ob.f_remove_olf_files();


  ob.f_temperatureCycle(mesh, mesh_neighbour, mesh_system_info, mesh_statistic_info, mesh_options, mesh_options.b_temperatureCycle);

  //ob.f_routine(mesh, mesh_neighbour, mesh_system_info, mesh_statistic_info, mesh_options, 3, 0);
  //ob.f_temperatureCycle (mesh, mesh_neighbour, mesh_system_info, mesh_statistic_info, mesh_options, 0.1, 10.0, 0.2);
  //ob.f_metropolisMethod (mesh, mesh_neighbour, mesh_system_info, 3);
  //ob.f_temperatureCycle(mesh, mesh_neighbour, mesh_system_info, mesh_statistic_info, mesh_options, mesh_options.b_temperatureCycle);
  //ob.f_metropolisMethod(mesh, mesh_neighbour, mesh_system_info, mesh_statistic_info, mesh_options, mesh_options.metropolisMethod);

  ob.f_about_system_inf(mesh_system_info,2);

  //read_inputfile("test_file.dat", mesh_system_info, mesh_options);

  return 0;
}

int read_inputfile(std::string filename,
                   framework::Fs_general_info & system_info,
                   framework::s_options & system_option)
{

  std::ifstream inputfileflux(filename, std::ios::in);
  if (!inputfileflux.is_open()) {
    std::cerr << "Error abriendo \" " << filename << " \" para escribir en el." << std::endl;
    abort ();
  }

  std::map <int, std::string> map_test;                                         // empty map container
  map_test.insert(std::pair <int, std::string> (1, "SYSTEM_DIMENSION"));        // insert elements in random order
  map_test.insert(std::pair <int, std::string> (2, "AMOUN_OF_CYCLES"));         // insert elements in random order
  map_test.insert(std::pair <int, std::string> (3, "ISOTHER_TEMPERATURE"));     // insert elements in random order

  map_test.insert(std::pair <int, std::string> (4, "R_HS"));                    // insert elements in random order
  map_test.insert(std::pair <int, std::string> (5, "R_LS"));                    // insert elements in random order
  map_test.insert(std::pair <int, std::string> (6, "EXCHANGE_J"));              // insert elements in random order
  map_test.insert(std::pair <int, std::string> (7, "SPRIN_K"));                 // insert elements in random order
  map_test.insert(std::pair <int, std::string> (8, "DELTA_ENERGY_D"));          // insert elements in random order

  map_test.insert(std::pair <int, std::string> (9, "STAR_TEMPERATURE"));        // insert elements in random order
  map_test.insert(std::pair <int, std::string> (10, "RATE_TEMPERATURE"));       // insert elements in random order
  map_test.insert(std::pair <int, std::string> (11, "FINAL_TEMPERATURE"));      // insert elements in random order

  map_test.insert(std::pair <int, std::string> (12, "TEMPERATURE_CASE"));       // insert elements in random order
  map_test.insert(std::pair <int, std::string> (13, "ROUTINE_CASE"));           // insert elements in random order
  map_test.insert(std::pair <int, std::string> (14, "METROPOLIS_CASE"));        // insert elements in random order
  map_test.insert(std::pair <int, std::string> (15, "SWITCH_ENERGY_CASE"));     // insert elements in random order
  map_test.insert(std::pair <int, std::string> (16, "ISIN_MODEL"));             // insert elements in random order

  map_test.insert(std::pair <int, std::string> (17, "SCREEN_PRINT_RATE"));      // insert elements in random order
  map_test.insert(std::pair <int, std::string> (18, "FILE_PRINT_RATE"));        // insert elements in random order


  std::map <int, std::string> ::iterator itr;                     // printing map gquiz1
  /*!
    1"SYSTEM_DIMENSION"
    2"AMOUN_OF_CYCLES"
    3"ISOTHER_TEMPERATURE"

    4"R_HS"
    5"R_LS"
    6"EXCHANGE_J"
    7"SPRIN_K"
    8"DELTA_ENERGY_D"

    9"STAR_TEMPERATURE"
    10"RATE_TEMPERATURE"
    11"FINAL_TEMPERATURE"

    12"TEMPERATURE_CASE"
    13"ROUTINE_CASE"
    14"METROPOLIS_CASE"
    15"SWITCH_ENERGY"
    16"ISIN_MODEL"

    17"SCREEN_PRINT_RATE"
    18"FILE_PRINT_RATE"
  */


  double n_values_number = 0.0;
  int line_counter = 0;

  while (! inputfileflux.eof() ){

    line_counter++;

    // read in whole line
    std::string line;
    getline(inputfileflux,line);
    //std::cout << line.c_str() << std::endl;
    //std::cout<<line<<std::endl;
    // ignore blank lines
    std::string empty = "";
    if(line == empty) continue;

    std::string hash = "#%";
    std::string equal = "=";
    std::string space = " ";
    std::string command, for_remove;
    std::string n_values;

    if(line.at(0) != hash.at(0) && line.at(0) != hash.at(1) ) break;


    //std::cout<<line<<std::endl;
    //std::istringstream iss(line,std::istringstream::in);
    //std::stringstream iss(line);
    //double values;


    //if(line.at(0) != equal.at(0) && line.at(0) != space.at(0) && !isdigit(line.at(0))){
    //    iss >> command >> for_remove >> values;
    //    std::cout << command<<" "<< for_remove<<" "<< values<<std::endl;
    //  }




     /*
      if(line.at(0) != hash.at(0) && line.at(0) != hash.at(0) && line.at(0) != equal.at(0) && line.at(0) != space.at(0) && !isdigit(line.at(0))){
              std::cout<<line<<std::endl;}
   */

 /*
    int i = 0;
    for( ; line.at(i) != equal.at(0) && line.at(i) != space.at(0) ;i++){

        if(line.at(i) == hash.at(0)) break;
        command.push_back(line.at(i));

      //else if (isdigit(line.at(i))) n_values.push_back(line.at(i));
    }
    for (;i < line.length(); ++i){
      if (isdigit(line.at(i))) {
        n_values.push_back(line.at(i));

        break;
      }

    }
    //std::cout<<"Numero"<< n_values<<std::endl;
    //std::string::size_type sz;     // alias of size_t
    //if (isdigit(n_values[0]) //n_values_number = std::stod(n_values, &sz);
    */
    int opt = 0;

    for (itr = map_test.begin(); itr != map_test.end(); ++itr) {
      if ( itr->second == command)
      {
        opt = itr ->first;
        //std::cout  <<  '\t' << itr->first
        //          <<  '\t' << itr->second << '\n';
      }

    }

    //std::cout<<"command:"<< command << std::endl;
    //std::cout<<"number:"<< n_values << std::endl;

 /*struct s_options{
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

    */
    switch(opt){

      case 1:{ //!SYSTEM_DIMENSION
        //system_info.FsNi =
        break;
      }

      case 2:{ //! AMOUN_OF_CYCLES
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 3:{//!ISOTHER_TEMPERATURE
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 4:{//!R_HS
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 5:{//!R_LS
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 6:{ //!EXCHANGE_J
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 7:{//!SPRIN_K
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 8:{//!DELTA_ENERGY_D
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 9:{//!STAR_TEMPERATURE
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 10:{//!RATE_TEMPERATURE
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 11:{ //!FINAL_TEMPERATURE
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 12:{//!TEMPERATURE_CASE
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 13:{//!ROUTINE_CASE
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 14:{//!METROPOLIS_CASE
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 15:{//!SWITCH_ENERGY
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 16:{ //!ISIN_MODEL
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 17:{//!SCREEN_PRINT_RATE
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      case 18:{//!FILE_PRINT_RATE
        //std::cout<<"estoy en:"<< opt << std::endl;
        break;
      }

      //default: std::cout<<"a esto a´un le falta"<<std::endl;
    }
  }




  inputfileflux.close();
  return 0;
}
