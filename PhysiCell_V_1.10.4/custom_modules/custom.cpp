/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
    Phenotype& phenotype = pC -> phenotype; 
    if (parameters.bools("enable_file_loading")){
    
    // Read the xml document and make sure document is correctly read;
    pugi::xml_document checkpointing_document;
    pugi::xml_node physicell_checkpoint_root;
    std::string filename = "./output_1/output00000100.xml";
    checkpointing_document.load_file(filename.c_str());
        
    pugi::xml_parse_result result = checkpointing_document.load_file( filename.c_str());
    if( result.status != pugi::xml_parse_status::status_ok )
	{
		std::cout << "Error loading " << filename << "!" << std::endl; 
		
	} else {
        std::cout << "Succesfully loaded " << filename << "!" << std::endl; 
    }
    
    // Get the root of file    
    physicell_checkpoint_root = checkpointing_document.child("MultiCellDS");
    std::cout << "And the root output is this"<<std::endl; 
    std::cout << physicell_checkpoint_root<<std::endl;
    // Get name of the root node 
    std::string name = xml_get_my_name( physicell_checkpoint_root );
    std::cout << "And the root name"<<std::endl; 
    std::cout << name<<std::endl;
    // Now see if I can find the labels node; 
    
    pugi::xml_node cellular_information_node = physicell_checkpoint_root.child("cellular_information");
    pugi::xml_node cell_populations_node = cellular_information_node.child("cell_populations");
    pugi::xml_node cell_population_node = cell_populations_node.child("cell_population");
    pugi::xml_node custom_node = cell_population_node.child("custom");
    pugi::xml_node simplified_data_node = custom_node.child("simplified_data");
    pugi::xml_node labels_node = simplified_data_node.child("labels"); 
    std::string name2 = xml_get_my_name( labels_node );
    std::cout << "And the labels name"<<std::endl; 
    std::cout << name2 <<std::endl;
    
    // Now loop through all labels node and see if I can find all the paramet names in the xml file
    pugi::xml_node node1 = labels_node.first_child(); 
	/* 
    int i = 0; 
	while( node1 )
	{
		std::string name = xml_get_my_name( node1 );
		std::string value = xml_get_my_string_value( node1 );
        std::string index = node1.attribute("index").value();
        std::string size = node1.attribute("size").value(); 
        
        std::cout<<value<<std::endl;
        std::cout<<index<<std::endl;           
        std::cout<<size<<std::endl;
                        
		node1 = node1.next_sibling(); 
                
		i++; 
	}
	*/
	
    // Read matlab cell data into a matrix 
    std::string matlab_filename = "./output_1/output00000100_cells.mat";
    std::vector<std::vector<double>> B = read_matlab( matlab_filename );
    
    // Get number of cells in the simulation
    int number_of_cells = B[0].size();
    
    // Loop through every cell and initiate it in the simulation
    for (int i = 0; i < number_of_cells; i++) { 
        // get cell definition and create a cell
        Cell_Definition* pCD = cell_definitions_by_index[B[6][i]]; 
        pC = create_cell( *pCD ); 
        // assig cell position
        node1 = labels_node.first_child(); 
            int j = 0; 
            while( node1 )
            {
                std::string value = xml_get_my_string_value( node1 );
                std::string index_str = node1.attribute("index").value();
                int index = std::stoi(index_str);
                std::string size = node1.attribute("size").value(); 
                
                if (value == "ID") {
                    pC->ID = B[index][i];
                } else if (value == "parent_ID") {
                    pC->parent_ID = B[index][i];
                } else if (value == "position") {
                    pC->assign_position({B[index][i],B[index+1][i],B[index+2][i]});
                } else if (value == "total_volume") {
                    pC->set_total_volume(B[index][i]);
                } else if (value == "current_phase") {
                    /*
                    std::cout<<"code"<<std::endl;
                    std::cout<<pC->phenotype.cycle.current_phase().code<<std::endl;
                    std::cout<<"B:"<<B[index][i]<<std::endl;
                    while (pC->phenotype.cycle.current_phase().code != B[index][i]){
                        // check which cycle model this is and how many phases are there.
                        // change the transition rates to 9e99, advance the cycle and rewrite the real transition rates at the end from matlab
                        pC->phenotype.cycle.data.transition_rate(0,1)=9e99;
                        pC->phenotype.cycle.data.transition_rate(1,2)=9e99;
                        pC->phenotype.cycle.data.transition_rate(2,3)=9e99;
                        pC->phenotype.cycle.data.transition_rate(3,0)=9e99;
                        pC->phenotype.cycle.advance_cycle(pC,pC->phenotype, 0.01);
                        
                        std::cout<<"New code"<<std::endl;
                        std::cout<<pC->phenotype.cycle.current_phase().code<<std::endl;
                    }*/
                } else if (value == "elapsed_time_in_phase") {
                    pC->phenotype.cycle.data.elapsed_time_in_phase = B[index][i]; 
                } else if (value == "nuclear_volume") {
                    pC->phenotype.volume.nuclear = B[index][i];
                } else if (value == "cytoplasmic_volume") {
                    pC->phenotype.volume.cytoplasmic = B[index][i];
                } else if (value == "fluid_fraction") {
                    pC->phenotype.volume.fluid_fraction = B[index][i];
                } else if (value == "calcification_fraction") {
                    pC->phenotype.volume.calcified_fraction = B[index][i];
                } else if (value == "orientation") {
                    pC->state.orientation[0] = B[index][i];
                    pC->state.orientation[1] = B[index+1][i];
                    pC->state.orientation[2] = B[index+2][i];
                } else if (value == "polarity") {
                    pC->phenotype.geometry.polarity = B[index][i];
                } else if (value == "velocity") {
                    //pC->set_velocity({B[index][i];B[index+1][i];B[index+2][i]};
                } else if (value == "pressure") {
                    pC->state.simple_pressure = B[index][i]; 
                } else if (value == "number_of_nuclei") {
                    pC->state.number_of_nuclei = B[index][i];
                } else if (value == "damage") {
                    pC->state.damage = B[index][i];
                } else if (value == "total_attack_time") {
                    pC->state.total_attack_time = B[index][i];
                } else if (value == "contact_with_basement_membrane"){
                    pC->state.contact_with_basement_membrane = B[index][i];
                } else if (value == "current_cycle_phase_exit_rate"){
                    int phase_index = pC->phenotype.cycle.data.current_phase_index; 
                    pC->phenotype.cycle.data.exit_rate(phase_index) = B[index][i];
                } else if (value == "dead") {
                    pC->phenotype.death.dead = B[index][i];
                } else if (value == "current_death_model"){
                    pC->phenotype.death.current_death_model_index = B[index][i];
                } else if (value == "death_rates"){
                    // Update apoptotic and necrotic death rates if these changed
                } else if (value == "cytoplasmic_biomass_change_rate"){
                    //pC->set_cytoplasmic_biomass_change_rate(B[index][i]);
                } else if (value == "nuclear_biomass_change_rate"){
                    //pC->set_nuclear_biomass_change_rate(B[index][i]);
                } else if (value == "fluid_change_rate"){
                    //pC->set_fluid_change_rate(B[index][i]);
                } else if (value == "calcification_rate"){
                    //pC->set_calcification_rate(B[index][i]);
                } else if (value == "target_solid_cytoplasmic"){
                    //pC->set_target_solid_cytoplasmic(B[index][i]);
                } else if (value == "target_solid_nuclear"){
                    //pC->set_target_solid_nuclear(B[index][i]);
                } else if (value == "target_fluid_fraction"){
                    //pC->set_target_fluid_fraction(B[index][i]);
                } else if (value == "radius"){
                    //pC->set_radius(B[index][i]);
                } else if (value == "nuclear_radius"){
                    //pC->set_nuclear_radius(B[index][i]);
                } else if (value == "surface_area"){
                    //pC->set_surface_area(B[index][i]);
                } else if (value == "cell_cell_adhesion_strength"){
                    //pC->set_cell_cell_adhesion_strength(B[index][i]);
                } else if (value == "cell_BM_adhesion_strength"){
                    //pC->set_cell_BM_adhesion_strength(B[index][i]);
                } else if (value == "cell_cell_repulsion_strength"){
                    //pC->set_cell_cell_repulsion_strength(B[index][i]);
                } else if (value == "cell_BM_repulsion_strength"){
                    //pC->set_cell_BM_repulsion_strength(B[index][i]);
                } else if (value == "cell_adhesion_affinities"){
                    //pC->set_cell_adhesion_affinities(B[index][i],size);
                } else if (value == "target_solid_nuclear"){
                    //pC->set_target_solid_nuclear(B[index][i]);
                } else if (value == "secretion_rates"){
                        pugi::xml_node variables_node = physicell_checkpoint_root.child("microenvironment").child("domain").child("variables");
                        pugi::xml_node substrate_node = variables_node.first_child();
                    for (int num_substrates = 0; num_substrates<std::stoi(size); num_substrates++){
                                                                        
                        std::string substrate_name = substrate_node.attribute("name").value(); 
                        double value = B[index+num_substrates][i];
                        set_single_behavior(pC,substrate_name+" secretion",value);
                            
                        substrate_node = substrate_node.next_sibling();
                       
                    }
                } else if (value == "transformation_rates"){
                    for (int cell_type_ID = 0; cell_type_ID<std::stoi(size); cell_type_ID++){
                        
                        double value = B[index+cell_type_ID][i];
                        set_single_behavior(pC,"transform to cell type "+std::to_string(cell_type_ID),value); 
                        
                    }
                } else {
                    
                    // Check if the value is in custom data
                    // Custom scalar
                    
                    for( int j=0 ; j < pC->custom_data.variables.size(); j++ )
                        {
                            std::string name = pC->custom_data.variables[j].name; 
                            if (value == name) {
                                pC->custom_data[value] = B[index+j][i]; 
                            }
                        }
                    // Custom vector 
                    for( int j=0 ; j < pC->custom_data.vector_variables.size(); j++ )
                        {
                            std::string name = pC->custom_data.vector_variables[j].name; 
                            if (value == name) {
                                for (int k = 0; k<std::stoi(size);k++){
                                pC->custom_data[value] = B[index+j+k][i]; 
                                }
                            }
                        }
                    
                }
                    
                    
                //std::cout<<value<<std::endl;
                //std::cout<<index<<std::endl;           
                //std::cout<<size<<std::endl;
                                
                node1 = node1.next_sibling(); 
                        
                j++; 
            }

    }
	
	/*
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	*/
    
    
        
       //load_from_checkpoint(parameters.strings("checkpoint_filename"));
        
    } else {
        
    Cell_Definition* pCD = find_cell_definition("wt"); 
    pC = create_cell( *pCD ); 
    pC->assign_position( {0.0,0.0,0.0} );
    
    }
    
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void load_from_checkpoint(std::string filename)
{
    // Read matlab cell data into a matrix 
    std::vector<std::vector<double>> B = read_matlab( filename );
    
    // Get number of cells in the simulation
    int number_of_cells = B[0].size();
    Cell* pC;
    // Loop through every cell and initiate it in the simulation
    for (int i = 0; i < number_of_cells; i++) { 
        // get cell definition and create a cell
        Cell_Definition* pCD = cell_definitions_by_index[B[6][i]]; 
        pC = create_cell( *pCD ); 
        // assig cell position
        pC -> assign_position({B[2][i],B[3][i],B[4][i]}); 
        
        // Get cell the same ID and parent ID
        pC -> ID = B[0][i];
        pC -> parent_ID = B[1][i];
        
        // Get cell volume
        pC -> set_total_volume( B[5][i] );
    }
    return;    
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 
