
#include "System.h"
#include "Storage.h"

Storage::Storage(
	std::weak_ptr<System> system_) {
	//std::cout << "FDM constructor" << std::endl;
	
	system = system_;

	std::ofstream statesOutput("Temp.sta");
	

	std::shared_ptr<System> SYSTEM = system.lock();//upgrades weak to shared

	
	if (SYSTEM) {
	 	statesOutput << "node_count " << SYSTEM->coordInfoVecs.nodeLocX.size() << '\n';
	 	statesOutput << "edge_count " << SYSTEM->coordInfoVecs.num_edges << '\n';
	 	statesOutput << "elem_count " << SYSTEM->coordInfoVecs.num_triangles << '\n';
	}


	statesOutput.close();
}


void Storage::print_VTK_File(void) {

	std::shared_ptr<System> SYSTEM = system.lock();

	if ((SYSTEM)) {

		iteration+=1;
		int digits = ceil(log10(iteration));// + 1));
		std::string format = ".vtk";
		std::string Number;
	  //std::string initial = "Animation_realistic/FINALLY_FLAT_";
    //std::string initial = "Animation_realistic/Simulations_with_flat_sheet_tests/some_boundary_tests_";
		//std::string initial = "Animation_realistic/Double_Sheet_Testing/10_10_24_double_tau_0_kbt_0_01505_lin_1_9_area_1_bend_15_LJep_0_";
    std::string initial = "Animation_realistic/Quals_sims/_Does_this_work_test_wl3-0hapf_";
    //1Vol_termination_off_Boundary_Updated_Weak_Region_Changed_Growth_Algorithm_on_8_26_24_changes_tracker_";
		// std::string initial = "Animation_realistic/YB_recalibedAlt_L0d75B0d135_test1_5thtry_";
		//std::string initial = "Regular_Simulation";
		// std::string initial = "Animation_realistic2/YB_recalibed_L0d75A0d75B0d135_GrowthFreq25_strain0d05_Hillcoef16_delay2_";
		
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration-1);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration-1);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration-1);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration-1);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());
		
	
		int numParticles = SYSTEM->generalParams.maxNodeCount;//one for lj particle

		
		

		
		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
		
		 
		ofs << "POINTS " <<numParticles + 1 << " float" << std::endl;
		for (int i = 0; i< numParticles; i++) {
			double xPos = SYSTEM->coordInfoVecs.nodeLocX[i];
			double yPos = SYSTEM->coordInfoVecs.nodeLocY[i];
			double zPos = SYSTEM->coordInfoVecs.nodeLocZ[i];
			
			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}
		ofs << std::setprecision(5) <<std::fixed<< SYSTEM->ljInfoVecs.LJ_PosX << " " << SYSTEM->ljInfoVecs.LJ_PosY << " " << SYSTEM->ljInfoVecs.LJ_PosZ << " " << '\n'<< std::fixed;
		
		

		
		int numEdges = SYSTEM->generalParams.true_num_edges;//coordInfoVecs.num_edges;//num_edges;
		int numCells = numEdges + 1;//one cell for LJ Particle, rest for edges of polymer
		int numNumsInCells = (3 * numEdges) + (2);//add one for lj and one to list it.
		
		
		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;
		//place edges as cells of type 2. 
		for (int edge = 0; edge < SYSTEM->coordInfoVecs.num_edges; edge++ ){
			if (SYSTEM->coordInfoVecs.edges2Nodes_1[edge] != INT_MAX &&
			 SYSTEM->coordInfoVecs.edges2Nodes_2[edge] != INT_MAX &&
			  SYSTEM->coordInfoVecs.edges2Nodes_1[edge] != -INT_MAX &&
			   SYSTEM->coordInfoVecs.edges2Nodes_2[edge] != -INT_MAX){
			int idA = SYSTEM->coordInfoVecs.edges2Nodes_1[edge];
			int idB = SYSTEM->coordInfoVecs.edges2Nodes_2[edge];
			

			ofs<< 2 << " " << idA << " " << idB << std::endl;
			}
		}

		ofs<< 1 << " " << SYSTEM->generalParams.maxNodeCount << std::endl;
		
		ofs << "CELL_TYPES " << numCells << std::endl;  
		//set edges and last set scattered points(dpd)
		for (int i = 0; i < (numEdges); i++) {
				
			ofs << 3 << std::endl;//edge (2d line)
		}
		ofs << 1 << std::endl;//edge (2d line)
		
	
		ofs << "CELL_DATA " << numCells << std::endl;
		ofs << "SCALARS Strain double " << std::endl;
		ofs << "LOOKUP_TABLE default "  << std::endl;
		//set strain for each edge

		for (int edge = 0; edge < SYSTEM->coordInfoVecs.num_edges; edge++ ){

			int idA = SYSTEM->coordInfoVecs.edges2Nodes_1[edge];
			int idB = SYSTEM->coordInfoVecs.edges2Nodes_2[edge];

			//if (idA >= SYSTEM->generalParams.maxNodeCount && idA != INT_MAX)
			//	std::cout<<idA<<std::endl;
			//if (idB >= SYSTEM->generalParams.maxNodeCount && idB != INT_MAX)
			//	std::cout<<idB<<std::endl;
			if (idA == INT_MAX || idB == INT_MAX || idA < 0 || idB < 0){
				continue;
			}
			double L0 = SYSTEM->generalParams.Rmin;
			double xL = SYSTEM->coordInfoVecs.nodeLocX[idA];
			double yL = SYSTEM->coordInfoVecs.nodeLocY[idA];
			double zL = SYSTEM->coordInfoVecs.nodeLocZ[idA];
			double xR = SYSTEM->coordInfoVecs.nodeLocX[idB];
			double yR = SYSTEM->coordInfoVecs.nodeLocY[idB];
			double zR = SYSTEM->coordInfoVecs.nodeLocZ[idB];
			
			double L1 = sqrt( (xL - xR)*(xL - xR)+(yL - yR)*(yL - yR)+(zL - zR)*(zL - zR));
			double strain = (L1 - L0) / L0;
			ofs << std::fixed << strain   << std::endl;
				
			
		}
		ofs << std::fixed << 0.1   << std::endl;
			
		ofs.close();
	
	}

	// if ((SYSTEM)) {

	// 	iteration_daughter+=1;
	// 	int digits = ceil(log10(iteration + 1));
	// 	std::string format = ".vtk";
	// 	std::string Number;
	// 	//std::string initial = "Animation_new/New_MC_interval_";
	// 	//std::string initial = "Animation_realistic_anneal/attempt1_kT3d0anneal_ks7d2_kb15_adh15_dt0d0002_";
	// 	//std::string initial = "Animation_QN/ss0d025_kT_0d5_QN_";
	// 	//std::string initial = "Animation_realistic_kT1d0/atp1_kT1d0_ks7d2_kb15_adh3d75_dt0d0002_";
	// 	//std::string initial = "Animation_5nm/atp1_kT2d5_ks28d8_kb15_ka720_adh15_dt0d0002_";
	// 	//std::string initial = "Animation_realistic_finaltry/wrap_v0d0001_dt0d0001_newrange_";
	// 	//std::string initial = "Animation_realistic/membrane_";//volumetest40_n2d0lowhem10_ka0_eqvol1d5_";//spheretest_rad0d17549_lowerhem5_ka5_ks25kb5_LJR2_"; //Anneal_adh15_Rv0d75_MD20a7d5_v0d2_NKBT4000_dt0d0002_";
	// 	//std::string initial = "Animation_realistic/yeastbudding_septinring_test_3particle_";
	// 	//std::string initial = "Animation_realistic/YB_cellwall3_l0d2b0d05a0d2_edgegrowth_areastrain0d2_maxt400_esfreq100_longsim200_";//yeastbudding_septin40_test_6particle_1pullonly_";
	// 	//std::string initial = "Animation_realistic9/YB_cellwall4_l1d0b0d1_st0d05_dampedring400_dt0d001_tau25_timecalib_damp1_ss0d05_formesh_";
	// 	std::string initial = "Animation_realistic9/YB_cellscission_daughter_";
	// 	//std::string initial = "Animation_realistic6/YB_cellwall3_l1d0b0d2a1d0_kT0d07_detgrowth_strain0d1_septin200_surfacegrowth_";
	// 	//std::string initial = "Animation_realistic/YB_reactive_isotropic_scalels0d01bs0d01as0d01_vols4d0_expthresh1d5rule1_";
	// 	//std::string initial = "Animation_realistic4/YB_cellwall3_isotropic_scalel0d5b0d25a0d5_vols0d3_expprob0d0025_dt0d001_1sttry_";
	// 	//std::string initial = "Animation_realistic/YB_recalib_l2a2b0d5_MPa0d0125125_";
	// 	//std::string initial = "Animation_realistic_flow/Pflow0d5_v0d0005_MRT0d005_dt0d0002_";
	// 	std::ofstream ofs;
	// 	if (digits == 1 || digits == 0) {
	// 		Number = "0000" + std::to_string(iteration_daughter);
	// 	}
	// 	else if (digits == 2) {
	// 		Number = "000" + std::to_string(iteration_daughter);
	// 	}
	// 	else if (digits == 3) {
	// 		Number = "00" + std::to_string(iteration_daughter);
	// 	}
	// 	else if (digits == 4) {
	// 		Number = "0" + std::to_string(iteration_daughter);
	// 	}

	// 	std::string Filename = initial + Number + format;

	// 	ofs.open(Filename.c_str());
		
	
	// 	int numParticles = SYSTEM->generalParams_daughter.maxNodeCount;//one for lj particle

		
		

		
	// 	ofs << "# vtk DataFile Version 3.0" << std::endl;
	// 	ofs << "Point representing Sub_cellular elem model" << std::endl;
	// 	ofs << "ASCII" << std::endl << std::endl;
	// 	ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
		
		 
	// 	ofs << "POINTS " <<numParticles + 1 << " float" << std::endl;
	// 	for (int i = 0; i< numParticles; i++) {
	// 		double xPos = SYSTEM->coordInfoVecs_daughter.nodeLocX[i];
	// 		double yPos = SYSTEM->coordInfoVecs_daughter.nodeLocY[i];
	// 		double zPos = SYSTEM->coordInfoVecs_daughter.nodeLocZ[i];
			
	// 		ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
	// 	}
	// 	ofs << std::setprecision(5) <<std::fixed<< SYSTEM->ljInfoVecs.LJ_PosX << " " << SYSTEM->ljInfoVecs.LJ_PosY << " " << SYSTEM->ljInfoVecs.LJ_PosZ << " " << '\n'<< std::fixed;
		
		

		
	// 	int numEdges = SYSTEM->generalParams_daughter.true_num_edges;//coordInfoVecs.num_edges;//num_edges;
	// 	int numCells = numEdges + 1;//one cell for LJ Particle, rest for edges of polymer
	// 	int numNumsInCells = (3 * numEdges) + (2);//add one for lj and one to list it.
		
		
	// 	ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;
	// 	//place edges as cells of type 2. 
	// 	int idA, idB;
	// 	//for (int edge = 0; edge < SYSTEM->coordInfoVecs_daughter.edges2Nodes_1.size(); edge++ ){
	// 	for (int edge = 0; edge < SYSTEM->coordInfoVecs_daughter.num_edges; edge++ ){
	// 		if (SYSTEM->coordInfoVecs_daughter.edges2Nodes_1[edge] < (INT_MAX-100) &&
	// 		 	SYSTEM->coordInfoVecs_daughter.edges2Nodes_2[edge] < (INT_MAX-100) &&
	// 		  	SYSTEM->coordInfoVecs_daughter.edges2Nodes_1[edge] >= 0 &&
	// 		   	SYSTEM->coordInfoVecs_daughter.edges2Nodes_2[edge] >= 0){
	// 			idA = SYSTEM->coordInfoVecs_daughter.edges2Nodes_1[edge];
	// 			idB = SYSTEM->coordInfoVecs_daughter.edges2Nodes_2[edge];}
	// 		else if (SYSTEM->coordInfoVecs_daughter.edges2Nodes_1[edge] < (INT_MAX-100) &&
	// 		  	SYSTEM->coordInfoVecs_daughter.edges2Nodes_1[edge] >= 0){
	// 			idA = SYSTEM->coordInfoVecs_daughter.edges2Nodes_1[edge];
	// 			idB = SYSTEM->coordInfoVecs_daughter.edges2Nodes_2[edge];}
	// 		else if (SYSTEM->coordInfoVecs_daughter.edges2Nodes_2[edge] < (INT_MAX-100) &&
	// 		   	SYSTEM->coordInfoVecs_daughter.edges2Nodes_2[edge] >= 0){
	// 			idA = SYSTEM->coordInfoVecs_daughter.edges2Nodes_1[edge];
	// 			idB = SYSTEM->coordInfoVecs_daughter.edges2Nodes_2[edge];}
	// 		else{continue;}

	// 		ofs<< 2 << " " << idA << " " << idB << std::endl;
	// 		//}
	// 	}

	// 	ofs<< 1 << " " << SYSTEM->generalParams_daughter.maxNodeCount << std::endl;
		
	// 	ofs << "CELL_TYPES " << numCells << std::endl;  
	// 	//set edges and last set scattered points(dpd)
	// 	for (int i = 0; i < (numEdges); i++) {
				
	// 		ofs << 3 << std::endl;//edge (2d line)
	// 	}
	// 	ofs << 1 << std::endl;//edge (2d line)
		
	
	// 	ofs << "CELL_DATA " << numCells << std::endl;
	// 	ofs << "SCALARS Strain double " << std::endl;
	// 	ofs << "LOOKUP_TABLE default "  << std::endl;
	// 	//set strain for each edge

	// 	for (int edge = 0; edge < SYSTEM->coordInfoVecs_daughter.num_edges; edge++ ){

	// 		int idA = SYSTEM->coordInfoVecs_daughter.edges2Nodes_1[edge];
	// 		int idB = SYSTEM->coordInfoVecs_daughter.edges2Nodes_2[edge];

	// 		//if (idA >= SYSTEM->generalParams.maxNodeCount && idA != INT_MAX)
	// 		//	std::cout<<idA<<std::endl;
	// 		//if (idB >= SYSTEM->generalParams.maxNodeCount && idB != INT_MAX)
	// 		//	std::cout<<idB<<std::endl;
	// 		if (idA >= (INT_MAX-100) || idB >= (INT_MAX-100) || idA < 0 || idB < 0){
	// 			continue;
	// 		}
	// 		double L0 = SYSTEM->generalParams.Rmin;
	// 		double xL = SYSTEM->coordInfoVecs_daughter.nodeLocX[idA];
	// 		double yL = SYSTEM->coordInfoVecs_daughter.nodeLocY[idA];
	// 		double zL = SYSTEM->coordInfoVecs_daughter.nodeLocZ[idA];
	// 		double xR = SYSTEM->coordInfoVecs_daughter.nodeLocX[idB];
	// 		double yR = SYSTEM->coordInfoVecs_daughter.nodeLocY[idB];
	// 		double zR = SYSTEM->coordInfoVecs_daughter.nodeLocZ[idB];
			
	// 		double L1 = sqrt( (xL - xR)*(xL - xR)+(yL - yR)*(yL - yR)+(zL - zR)*(zL - zR));
	// 		double strain = (L1 - L0) / L0;
	// 		ofs << std::fixed << strain   << std::endl;
				
			
	// 	}
	// 	ofs << std::fixed << 0.1   << std::endl;
			
	// 	ofs.close();
	
	// }

	//now print out the file for the capsid
	/*if ((SYSTEM)) {
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		//std::string initial = "Animation_realistic/yeastbudding_septinring_nucleus_test_3particle_";
		std::string initial = "Animation_realistic6/withseptin_rule3_lw15bw7d5aw15_volsp15_particle_dt0d00002_";//yeastbudding_septin40_nucleus_test_6particle_1pullonly_";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());
		
	
		unsigned numParticles = SYSTEM->generalParams.maxNodeCountLJ;

		unsigned num_connections=0;
		num_connections = SYSTEM->generalParams.maxNodeCountLJ;
		
		
		double xPos;
		double yPos;
		double zPos;
		
		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
		
		 
		ofs << "POINTS " <<numParticles + num_connections << " float" << std::endl;
		for (unsigned i = 0; i< numParticles; i++) {
			xPos = SYSTEM->ljInfoVecs.LJ_PosX_all[i];
			yPos = SYSTEM->ljInfoVecs.LJ_PosY_all[i];
			zPos = SYSTEM->ljInfoVecs.LJ_PosZ_all[i];
			
			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}
		//std::cout<<'here'<<std::flush;

		//set location for nodes that capside is connected to
		//ie  
		for (unsigned i = 0; i < num_connections; i++ ) {
			unsigned mem_id = i;
			//std::cout<<" "<< std::endl;
			//std::cout<<mem_id<<std::flush;
			xPos = SYSTEM->ljInfoVecs.LJ_PosX_all[mem_id];
			yPos = SYSTEM->ljInfoVecs.LJ_PosY_all[mem_id];
			zPos = SYSTEM->ljInfoVecs.LJ_PosZ_all[mem_id];
			
			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		
		}


		//std::cout<<'here1'<<std::flush;
		
		unsigned numCells = 1;
		numCells += num_connections;//add conections cells for edges

		unsigned numNumsInCells = 1 + numParticles;
		numNumsInCells += 3 * num_connections;//3 numbers per edge

		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;
		//place edges as cells of type 2. 
		ofs<< numParticles << " ";
		for (unsigned point = 0; point < numParticles; point++ ){
			ofs<< " " << point;
		}
		ofs<<" "<< std::endl;

		//std::cout<<'here2'<<std::flush;
		for (unsigned edge = 0; edge < num_connections; edge++ ){

			unsigned mem_id = edge;//numParticles + edge;
			unsigned cap_id = edge;//SYSTEM->capsidInfoVecs.tempCapsideId[edge];
				
			ofs <<2<< " "<< mem_id << " "<< cap_id <<std::endl;
		}
		ofs << "CELL_TYPES " << numCells << std::endl;  
		//set edges and last set scattered points
				
		ofs << 2 << std::endl;//scatter points for capsid
		
	//	std::cout<<'here3'<<std::flush;
		for (unsigned edge = 0; edge< num_connections; edge++ ){
			ofs<< 3 <<std::endl;
		}
		ofs.close();
	}*/
	
};

void Storage::storeVariables(void) {
	std::shared_ptr<System> SYSTEM = system.lock();
	if (SYSTEM) {

		iteration2+=1;
		int digits = ceil(log10(iteration2));// + 1));
		std::string format = ".sta";
		std::string Number;
		//std::string initial = "Animation_new/New_MC_interval_";
		//std::string initial = "Animation_realistic_anneal/attempt1_kT3d0anneal_ks7d2_kb15_adh15_dt0d0002_";
		//std::string initial = "Animation_QN/ss0d025_kT_0d5_QN_";
		//std::string initial = "Animation_realistic_kT1d0/atp1_kT1d0_ks7d2_kb15_adh3d75_dt0d0002_";
		//std::string initial = "Animation_5nm/atp1_kT2d5_ks28d8_kb15_ka720_adh15_dt0d0002_";
		//std::string initial = "Animation_realistic_finaltry/wrap_v0d0001_dt0d0001_newrange_";
		//std::string initial = "Animation_realistic/membrane_";//volumetest40_n2d0lowhem10_ka0_eqvol1d5_";//spheretest_rad0d17549_lowerhem5_ka5_ks25kb5_LJR2_"; //Anneal_adh15_Rv0d75_MD20a7d5_v0d2_NKBT4000_dt0d0002_";
		//std::string initial = "Animation_realistic/yeastbudding_septinring_test_3particle_";
		//std::string initial = "Variables_realistic/YB_reactive_isotropic_scalels0d01bs0d01as0d01_vols4d0_expthresh1d5rule1_";//yeastbudding_septin40_test_6particle_1pullonly_";
		//std::string initial = "Variables_realistic/YB_cellwall3_l0d2b0d05a0d2_edgegrowth_areastrain0d2_maxt400_esfreq100_longsim200_";
		std::string initial = "Variables_realistic/YB_cellwall4_newinitialmesh6_";
		//std::string initial = "Variables_realistic6/YB_cellwall3_l1d0b0d02a1d0_kT0d07_detgrowth_strain0d1_septin200_surfacegrowth_";
		//std::string initial = "Animation_realistic_flow/Pflow0d5_v0d0005_MRT0d005_dt0d0002_";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration2-1);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration2-1);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration2-1);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration2-1);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());
		//first create a new file using the current network strain
		
		
		/*for (int i = 0; i < SYSTEM->coordInfoVecs.num_edges; i++){//SYSTEM->coordInfoVecs.nodeLocX.size(); i++) {
			double x = SYSTEM->generalParams.angle_per_edge[i];
			
			ofs << std::setprecision(5) <<std::fixed<< x <<std::endl;
		
		}*/



		// double total_energy =  SYSTEM->linearSpringInfoVecs.linear_spring_energy + 
        // 						SYSTEM->areaTriangleInfoVecs.area_triangle_energy + 
        // 						SYSTEM->bendingTriangleInfoVecs.bending_triangle_energy + 
		// 						0.5*SYSTEM->linearSpringInfoVecs.memrepulsion_energy +
        // 						SYSTEM->ljInfoVecs.lj_energy_M +
		// 						SYSTEM->ljInfoVecs.lj_energy_LJ;

								
								
		//ofs << std::setprecision(5) <<std::fixed<< "total_energy=" << total_energy<<std::endl;
		//for (int i = 0; i < SYSTEM->ljInfoVecs.LJ_PosX_all.size(); i++){
		//	ofs << std::setprecision(5) <<std::fixed<< "nucleus " << SYSTEM->ljInfoVecs.LJ_PosX_all[i]<<" "<< SYSTEM->ljInfoVecs.LJ_PosY_all[i]<<" "<< SYSTEM->ljInfoVecs.LJ_PosZ_all[i]<<std::endl;
		//}
		unsigned numParticles = SYSTEM->generalParams.maxNodeCountLJ;
		
		ofs << std::setprecision(5) <<std::fixed<<"number of LJ particles "<<numParticles<<std::endl;
		
		// for (int i = 0; i < SYSTEM->coordInfoVecs.num_edges; i++){
		// 	ofs << std::setprecision(5) <<std::fixed<<"Angle = "<<SYSTEM->generalParams.angle_per_edge[i]<<std::endl;
		// }

		//place nodes
		for (int i = 0; i < SYSTEM->generalParams.maxNodeCount; i++){//SYSTEM->coordInfoVecs.nodeLocX.size(); i++) {
			double x = SYSTEM->coordInfoVecs.nodeLocX[i];
			double y = SYSTEM->coordInfoVecs.nodeLocY[i];
			double z = SYSTEM->coordInfoVecs.nodeLocZ[i];
			double nodes_in_tip = SYSTEM->generalParams.nodes_in_tip[i];
			double nodes_in_upperhem = SYSTEM->generalParams.nodes_in_upperhem[i];
			ofs << std::setprecision(5) <<std::fixed<< "<node>" << x << " " << y << " " << z <<"</node>"<<std::endl;//" "<<nodes_in_tip<<" "<<nodes_in_upperhem<<std::endl;
		
		}

		for (int i = 0; i < SYSTEM->coordInfoVecs.num_triangles; i++){//SYSTEM->coordInfoVecs.triangles2Nodes_1.size(); i++) {
			int t2n_1 = SYSTEM->coordInfoVecs.triangles2Nodes_1[i];
			int t2n_2 = SYSTEM->coordInfoVecs.triangles2Nodes_2[i];
			int t2n_3 = SYSTEM->coordInfoVecs.triangles2Nodes_3[i];
			ofs << std::setprecision(5) <<std::fixed<< "<elem> " << t2n_1+1 << " " << t2n_2+1 << " " << t2n_3+1 <<" </elem>"<<std::endl;
		
		}

		// //for (int i = 0; i < SYSTEM->generalParams.maxNodeCount; i++){
		// //	double t2n_1 = SYSTEM->coordInfoVecs.nodeForceX[i];
		// //	double t2n_2 = SYSTEM->coordInfoVecs.nodeForceY[i];
		// //	double t2n_3 = SYSTEM->coordInfoVecs.nodeForceZ[i];
		// //	ofs << std::setprecision(5) <<std::fixed<< " " << t2n_1 << " " << t2n_2 << " " << t2n_3 <<" "<<std::endl;
		// //
		// //}

		for (int i = 0; i < SYSTEM->coordInfoVecs.num_triangles; i++) {
			int t2n_1 = SYSTEM->coordInfoVecs.triangles2Edges_1[i];
			int t2n_2 = SYSTEM->coordInfoVecs.triangles2Edges_2[i];
			int t2n_3 = SYSTEM->coordInfoVecs.triangles2Edges_3[i];
			ofs << std::setprecision(5) <<std::fixed<< "<elem2edge> " << t2n_1+1 << " " << t2n_2+1 << " " << t2n_3+1 <<" </elem2edge>"<<std::endl;
		
		}

		for (int i = 0; i < SYSTEM->coordInfoVecs.num_edges; i++) {
			int t2n_1 = SYSTEM->coordInfoVecs.edges2Nodes_1[i];
			int t2n_2 = SYSTEM->coordInfoVecs.edges2Nodes_2[i];
			ofs << std::setprecision(5) <<std::fixed<< "<edgeinfo> " << t2n_1+1 << " " << t2n_2+1 <<" </edgeinfo>"<<std::endl;
		
		}

		for (int i = 0; i < SYSTEM->coordInfoVecs.num_edges; i++) {
			int t2n_1 = SYSTEM->coordInfoVecs.edges2Triangles_1[i];
			int t2n_2 = SYSTEM->coordInfoVecs.edges2Triangles_2[i];
			ofs << std::setprecision(5) <<std::fixed<< "<edge2elem> " << t2n_1+1 << " " << t2n_2+1 <<" </edge2elem>"<<std::endl;
		
		}



		for (int i = 0; i < SYSTEM->generalParams.maxNodeCount; i++) {
			int nn1 = SYSTEM->coordInfoVecs.nndata1[i];
			int nn2 = SYSTEM->coordInfoVecs.nndata2[i];
			int nn3 = SYSTEM->coordInfoVecs.nndata3[i];
			int nn4 = SYSTEM->coordInfoVecs.nndata4[i];
			int nn5 = SYSTEM->coordInfoVecs.nndata5[i];
			int nn6 = SYSTEM->coordInfoVecs.nndata6[i];
			int nn7 = SYSTEM->coordInfoVecs.nndata7[i];
			int nn8 = SYSTEM->coordInfoVecs.nndata8[i];
			int nn9 = SYSTEM->coordInfoVecs.nndata9[i];
			// int nn10 = SYSTEM->coordInfoVecs.nndata10[i];
			// int nn11 = SYSTEM->coordInfoVecs.nndata11[i];
			// int nn12 = SYSTEM->coordInfoVecs.nndata12[i];
			ofs << std::setprecision(5) <<std::fixed<< "<nndata> " << nn1+1 << " " << nn2+1 <<" "<< nn3+1 <<" "<< nn4+1 <<" "<< nn5+1 <<" "<< nn6+1 <<" "<< nn7+1 <<" "<< nn8+1 <<" "<< nn9+1 <<" </nndata> "<<std::endl;
		
		}

	}
}
