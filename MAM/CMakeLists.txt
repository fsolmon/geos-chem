add_library(MAM STATIC EXCLUDE_FROM_ALL
        mam_utils.F90
	physconst.F90
	chem_mods.F90
        constituents.F90
        phys_control.F90 
	modal_aero_data.F90
	modal_aero_amicphys.F90
        modal_aero_newnuc.F90
        modal_aero_coag.F90
	modal_aero_rename.F90
        modal_aero_gasaerexch.F90
	modal_aero_wateruptake.F90
	modal_aero_calcsize.F90
	modal_aero_initialize_data.F90
        module_data_mosaic_aero.F90
    	module_mosaic_box_aerchem.F90
      	module_mosaic_lsode.F90
        module_data_mosaic_asecthp.F90
      	module_mosaic_cam_init.F90
   	module_mosaic_support.F90
        module_data_mosaic_main.F90
    	module_mosaic_ext.F90
        module_mosaic_astem.F90
	module_mosaic_init_aerpar.F90
	radconstants.F90
        rad_constituents.F90 
	physics_buffer.F90
	physics_types.F90
        buffer.F90	
        wv_saturation.F90
	error_function.F90
	)
target_link_libraries(MAM
	PUBLIC GeosUtil
)

