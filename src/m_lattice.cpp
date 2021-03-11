#include<m_lattice.h>
#include<framework.h>


m_lattice::m_lattice(){}

m_lattice::~m_lattice(){}



int m_lattice::f_create_2D_lattice(std::vector<s_two_D_system> & system_in,
                                   m_lattice::s_general_info & system_info)
{
  std::cout <<"\n----------Creating 2D lattice with (poner aqui la variable del total de atomos) atoms.----------"<<std::endl;
  std::cout << "system_info.Fsdelta:" <<system_info.Fsdelta <<std::endl;
  
  /*!---------------- creating system and neighbour of external variable ----------------*/

  system_in.resize(system_info.FsTotalAtoms);           /*! creating the external vector system.*/
 
  
  /*!---------------- initialising system system and neighbour ----------------*/
  //double x_b = double (0.9 - (system_info.FsR_ls + system_info.FsR_ls));
  //f_det_VandL(system_in, system_info, 1);
 
  int flag = 0;
  for (int y = 0; y < system_info.FsNi; ++y){
    for (int x = 0; x < system_info.FsNi; ++x){
      system_in[flag].tSpin = 0;
      //system_in[flag].sBorder.resize(6); //int_randomNumber_generator(0, 1);
      //if (system_in[flag].tSpin == 1) system_in[flag].R = user_R_hs_d; //cambiar esto por una operacion matematica
      //else  system_in[flag].R = user_R_hs_d/1.1;
      system_in[flag].spin = 1;
      system_in[flag].tSpin = 0;
      system_in[flag].x = x * 0.5;
      system_in[flag].y = y * 0.5;
      system_in[flag].two_D_neighbour.resize(4);
      
      flag++;
    }
  }
  
  //std::cout << "system_info.Fsdelta:" <<system_info.Fsdelta <<std::endl;
  //system_info.FsL = f_det_VandL (system_in, system_info, 2);
  //system_info.FsL = system_info.Fsdelta * system_info.FsNi;
  //std::cout<< "system_info.FsL:"<<system_info.FsL<<std::endl;
  //std::cout << "system_info.FsL\tsystem_info.Fsdelta\tmain_system_info.FsNi\n" <<system_info.FsL <<"\t"<<system_info.Fsdelta <<"\t" <<main_system_info.FsNi<<std::endl;
  //system_info.FsVolume = system_info.FsL * system_info.FsL * system_info.FsL;
  std::puts("Done!");
  return 0;
}

int m_lattice::f_create_1D_lattice(std::vector<s_one_D_system> & system_in_one,
                                   s_general_info & system_info)
{
  std::cout <<"\n----------Creating 1D lattice with (poner aqui la variable del total de atomos) atoms.----------"<<std::endl;
  std::cout << "system_info_one.Fsdelta:" <<system_info.Fsdelta <<std::endl;
  
  /*!---------------- creating system and neighbour of external variable ----------------*/

  system_in_one.resize(system_info.FsTotalAtoms);           /*! creating the external vector system.*/
 
  for (int x = 0; x < system_info.FsNi; ++x){  /*! initializing the 1 dimension array */
    system_in_one[x].tSpin = 0;
    system_in_one[x].spin = 1;
    system_in_one[x].x = x * 0.5;
  }

  std::puts("Done!");
  return 0;
}

int m_lattice::run(std::vector<s_one_D_system> & system_in_one,
                   s_general_info & system_info)
{
  f_create_1D_lattice(system_in_one, system_info);
  return 0;
}