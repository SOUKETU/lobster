
# Load libraries
library(TMB)               
library(SpatialDeltaGLMM)

file.sources = list.files(path="/Users/jiecao/Desktop/Google Drive/work/Postdoc research/GOM_lobster_spatial_distribution/spatio-temporal-model-lobster/R",pattern="*.R",full.names=TRUE)
for (f in file.sources) {
  source(f)
  print(f)}

###############
# Settings
###############
covariate=0
sce=c("total","female","male","juvenile","adult")[1]

Data_Set = c("GOM_lobster")
Sim_Settings = list("Species_Set"=1:100, "Nyears"=10, "Nsamp_per_year"=600, "Depth_km"=-1, "Depth_km2"=-1, "Dist_sqrtkm"=0, "SigmaO1"=0.5, "SigmaO2"=0.5, "SigmaE1"=0.5, "SigmaE2"=0.5, "SigmaVY1"=0.05, "Sigma_VY2"=0.05, "Range1"=1000, "Range2"=500, "SigmaM"=1)
Version = "geo_index_v4b"
Method = c("Grid", "Mesh")[2]
grid_size_km = 25
n_x = c(100, 250, 500, 1000, 2000)[1] # Number of stations
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=1) # 1=Presence-absence; 2=Density given presence; #Epsilon=Spatio-temporal; #Omega=Spatial
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) # Structure for beta or epsilon over time: 0=None (default); 1=WhiteNoise; 2=RandomWalk; 3=Constant
VesselConfig = c("Vessel"=0, "VesselYear"=0)
ObsModel = 2  # 0=normal (log-link); 1=lognormal; 2=gamma; 4=ZANB; 5=ZINB; 11=lognormal-mixture; 12=gamma-mixture
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid

# Determine region
Region = switch( Data_Set,"GOM_lobster"="West_GOM", "SAWC_jacopever"="South_Africa", "Sim"="California_current")
strata.limits = list('All_areas'=c('a','b','c','d','e'),'a'=c('a'),'b'=c('b'),'c'=c('c'),'d'=c('d'),'e'=c('e'))

# This is where all runs will be located
if(covariate==1){
  DateFile = paste('/Users/jiecao/Desktop/Google Drive/work/Postdoc research/GOM_lobster_spatial_distribution/spatio-temporal-model-lobster/model results','/',Sys.Date(),'-','covariates-',sce,'/',sep='')
}else{
  DateFile = paste('/Users/jiecao/Desktop/Google Drive/work/Postdoc research/GOM_lobster_spatial_distribution/spatio-temporal-model-lobster/model results','/',Sys.Date(),'-','no-covariates-',sce,'/',sep='')
}
dir.create(DateFile)

# Save options for future records
Record = ThorsonUtilities::bundlelist( c("Data_Set","Sim_Settings","strata.limits","Region","Version","Method","grid_size_km","n_x","FieldConfig","RhoConfig","VesselConfig","ObsModel","Kmeans_Config") )
save( Record, file=file.path(DateFile,"Record.RData"))
capture.output( Record, file=paste0(DateFile,"Record.txt"))

################
# Prepare data
################

setwd('/Users/jiecao/Desktop/Google Drive/work/Postdoc research/GOM_lobster_spatial_distribution/spatio-temporal-model-lobster/data')

# Read in lobster survey data
load("LobsterSurveyData.rda")

if (sce=="total")
  lobster_catch = LobsterSurveyData$Catch_n
if (sce=="female")
  lobster_catch = LobsterSurveyData$Female2
if (sce=="male")
  lobster_catch = LobsterSurveyData$Male2
if (sce=="juvenile")
  lobster_catch = LobsterSurveyData$Juvenile2
if (sce=="adult")
  lobster_catch = LobsterSurveyData$Adult2

Data_Geostat = data.frame("Catch_N"=lobster_catch, "Year"=LobsterSurveyData[,'Year'], 
                          "Vessel"="missing", "AreaSwept_km2"=LobsterSurveyData[,'Tow_length_nm']*1.852*10.50563*0.001, 
                          "Lat"=LobsterSurveyData[,'Lat'], "Lon"=LobsterSurveyData[,'Lon'],"Depth"=LobsterSurveyData[,'Depth'])

Year_Set = sort(unique(Data_Geostat[,'Year']))
Q_ik = as.matrix(LobsterSurveyData[,'DOY'],nrow=length(LobsterSurveyData[,'DOY']),ncol=1)
# Get extrapolation data

if( Region == "West_GOM" ){
  Extrapolation_List = Prepare_WGOM_Extrapolation_Data_Fn( strata.limits=strata.limits )
}

Data_Geostat = na.omit( Data_Geostat )

# Calculate spatial information for SPDE mesh, strata areas, and AR1 process
Spatial_List = Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile )
Data_Geostat = cbind( Data_Geostat, Spatial_List$loc_UTM, "knot_i"=Spatial_List$knot_i )

################
# Make and Run TMB model
# (THIS WILL BE SIMILAR FOR EVERY DATA SET) 
################

# Make TMB data list
TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "b_i"=Data_Geostat[,'Catch_N'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year'], "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method, "Options"=c(SD_site_density=0, SD_site_logdensity=0, Calculate_Range=1, Calculate_evenness=0, Calculate_effective_area=1) )

# Make TMB object
TmbList = Build_TMB_Fn("TmbData"=TmbData, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "VesselConfig"=VesselConfig, "loc_x"=Spatial_List$loc_x)
Obj = TmbList[["Obj"]]

# Run model
Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=FALSE )
Report = Obj$report()

# Save stuff
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))

################
# Make diagnostic plots
################

# Plot settings
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

# Plot Anisotropy
PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Report, TmbData=TmbData )

# Plot surface
MapDetails_List = MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
PlotResultsOnMap_Fn(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.0)
                                                                                                                         
# Plot index
PlotIndex_Fn( DirName=DateFile, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, strata_names=ls(strata.limits), use_biascorr=TRUE )

# Plot center of gravity
Plot_range_shifts(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), FileName_COG=paste0(DateFile,"center_of_gravity.png"))

# Vessel effects
Return = Vessel_Fn(TmbData=TmbData, Sdreport=Opt[["SD"]], FileName_VYplot=paste0(DateFile,"VY-effect.jpg"))

# Positive catch rate Q-Q plot
Q = QQ_Fn( TmbData=TmbData, Report=Report, FileName_PP=paste0(DateFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(DateFile,"Posterior_Predictive-Histogram.jpg"), FileName_QQ=paste0(DateFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(DateFile,"Q-Q_hist.jpg"))


  