include ../make.inc

SOURCES=WENO2D_module.f90 \
	Riemann_Module.f90 \
	WENO_module.f90 \
	EcapTable_module.f90 \
	FlameTable_module.f90 \
	LevelSet3_module.f90 \
	Helmeos_module.f90 \
	Turb_module.f90 \
	PPT_module.f90 \
	NuSpec_module.f90 \
	GW_module.f90 \
	Main.f90 \
	AlphaSplit_r.f90 \
	AlphaSplit_z.f90 \
	BackUp_RVE.f90 \
	BackupRefArray.f90 \
	Boundary1D.f90 \
	Boundary1D_FULL.f90 \
	Boundary2D.f90 \
	Boundary2D_X.f90 \
	BuildWENO.f90 \
	BuildHydro.f90 \
	CheckRho.f90 \
	CloseFile.f90 \
	EOSTABLE.f90 \
	FermiMo.f90 \
	FindAsh.f90 \
	Final.f90 \
	FindCentralDensity.f90 \
	FindDt.f90 \
	FindEnergy.f90 \
	FindFlameVel.f90 \
	FindGamma.f90 \
	FindMass.f90 \
	FindPressure.f90 \
	FindPotential.f90 \
	FromRVEToU.f90 \
	GetEpsilon.f90 \
	GetGrid.f90 \
	GetNuQdot.f90 \
	GetRho_EOSPToR.f90 \
	GetRho_EOSRToP.f90 \
	GetRho.f90 \
	GetRho_template.f90 \
	GetSponge.f90 \
	GetVel.f90 \
	GetYe.f90 \
	Highorder.f90 \
	Helmeos.f90 \
	Ini_Der.f90 \
	Initial_x2.f90 \
	Interpolation.f90 \
	OpenFile.f90 \
	Output_log.f90 \
	Output_profile.f90 \
	PotentialRelax.f90 \
	RungeKutta.f90 \
	ReportParameter.f90 \
	ReportCurrentData.f90 \
	Spatial.f90 \
	Update.f90 

OBJECTS=$(SOURCES:.f90=.o )

WENO2D: WENO2D_module.o $(OBJECTS)  
	$(F90) $(LDFLAGS) -o WENO2D $(OBJECTS) 

$(OBJECTS): %.o: %.f90 
	$(F90) $(F90FLAGS) -c $< -o $@

weno2d.mod: WENO2D_module.f90
	$(F90) -C WENO2D_module.f90

weno.mod: WENO_module.f90
	$(F90) -C WENO_module.f90

riemann.mod: Riemann_Module.f90
	$(F90) -C Riemann_Module.f90

ecap.mod: EcapTable_module.f90
	$(F90) -C EcapTable_module.f90

flametable.mod: flametable_module.f90
	$(F90) -C FlameTable_module.f90

turb.mod: Turb_module.f90
	$(F90) -C Turb_module.f90

ppt.mod: PPT_module.f90
	$(F90) -C PPT_module.f90

nuspec.mod: NuSpec_module.f90
	$(F90) -C NuSpec_module.f90

levelset.mod: LevelSet3_module.f90
	$(F90) -C LevelSet3_module.f90

gwave.mod: GW_module.f90
	$(F90) -C GW_module.f90

clean:
	rm -rf WENO2D
	rm -rf *.o 	
	rm -rf *.mod
	rm tmp.txt
cleanfile:
	rm -rf ./Outfile/Star_WENO_*
	rm -rf ./Outfile/Hydro/Star_WENO_*
	rm -rf ./Outfile/Flame/Star_WENO_*
	rm -rf ./Outfile/Isotope/Star_WENO_*
	rm -rf ./Outfile/Neutrino/Star_WENO_*
	rm -rf ./Outfile/Turbulence/Star_WENO_*
	rm -rf ./Outfile/GWave/Star_WENO_*
	rm -rf ./Tracer/Star_WENO_*
cleanpic:
	rm -rf ./Plot/Hydro/*.png
	rm -rf ./Plot/Flame/*.png
	rm -rf ./Plot/Isotope/*.png
	rm -rf ./Plot/Neutrino/*.png
	rm -rf ./Plot/Turbulence/*.png
	rm -rf ./Plot/GWave/*.png
	rm -rf ./Plot/Tracer/*.png