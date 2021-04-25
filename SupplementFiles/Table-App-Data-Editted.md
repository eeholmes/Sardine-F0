---
output:
  pdf_document: default
  html_document: default
---
Landings Data
=============

<p>
We used landings (in metric tons) of oil sardines in
Kerala State 1956-2015. The data were collected and
processed into a total landings estimate based on a stratified sampling
of landing sites along the southwest coast of India throughout the year. The
program is run by the Central Marine Fisheries Research Institute
(CMFRI) in Cochin, India. We obtained the data from reports published
by CMFRI; see references.
</p>
<h3>
References
</h3>
<p>
CMFRI reports were downloaded from the CMFRI Publication repository
<a href="http://www.cmfri.org.in">http://www.cmfri.org.in</a>.
</p>
<p>
1956-1968 Antony Raja BT
(1969). “Indian oil sardine.” <em>CMFRI Bulletin</em>, <b>16</b>, 1-142.
</p>
<p>
1968-1978 Pillai VN (1982). <em>Physical characteristics of the coastal
waters off the south-west coast of India with an attempt to study the
possible relationship with sardine, mackerel and anchovy fisheries</em>.
Thesis, University of Cochin.
</p>
<p>
1975-1984 Kerala Jacob T, Rajendran V, Pillai PKM, Andrews J, Satyavan
UK (1987). “An appraisal of the marine fisheries in Kerala.” Central
Marine Fisheries Research Institute. Report.
</p>
<p>
1975-1984 Kernataka Kurup KN, Nair GKK, Annam VP, Kant A, Beena MR,
Kambadkar L (1987). “An appraisal of the marine fisheries of Karnataka
and Goa.” Central Marine Fisheries Research Institute. Report.
</p>
<p>
1985-2015 Provided by CMFRI directly via a data request.
</p>
<hr />


SST Data
========

<p>
We used two primary reanalysis SST data sets: ioSST and ICOADS SST. We also compared to a third non-reanalysis data set. Reanalysis means that the data set uses multiple data sources and may use interpolation to fill in missing values. See Supplement 1 for comparisons of the data sources on a grid. 
</p>
<p>
**AVHRR Data** For the paper, we used the Daily Optimum Interpolation (OI), AVHRR Only, Version 2.1 data set (oiSST) by the Group for High Resolution Sea Surface Temperature (GHRSST). "AVHRR Only" is in contrast to the "AMSR+AVHRR Product (AMSR = Advanced Microwave Scanning Radiometer).

<a href="https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.nodc:GHRSST-AVHRR_OI-NCEI-L4-GLOB/html">https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.nodc:GHRSST-AVHRR_OI-NCEI-L4-GLOB/html</a>

and downloaded from 

<a href="https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg.html">https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg.html</a>

This data set is on a 0.25 degree and provides complete ocean temperature fields constructed by combining bias-adjusted observations from different platforms (AVHRR satellite instruments, ships, buoys) on a regular global grid, with gaps filled in by interpolation. 
</p>
<p>The ioSST data were compared to non-reanalysis AVHRR data.
For 1981 to 2003, we used the Pathfinder Version 5.2 (L3C) monthly day
and night product on a 0.0417 degree grid. These SST data use the
Advanced Very-High Resolution Radiometer (AVHRR) instrument on the
Pathfinder satellites and  were provided by the 
Group for High Resolution Sea Surface Temperature (GHRSST) and the US
National Oceanographic Data Center. This project was supported in part
by a grant from the NOAA Climate Data Record (CDR) Program for
satellites.
For 2004 to 2016, we used the NOAA CoastWatch sea surface temperature
(SST) products derived from NOAA’s Polar Operational Environmental
Satellites (POES). The SST estimates use the Advanced Very-High
Resolution Radiometer (AVHRR) instruments on the POES satellites and are
on a 0.1 degree grid.
</p>
<p>
**ICOADS** The International Comprehensive Ocean-Atmosphere Data Set (ICOADS) is a collection of surface marine data. SST data from 1960 onward were used, which are on a 1 degree grid. The nearshore data (boxes 1 to 5) are not as accurate as the AVHRR-based SST data. The ICOADS data were only used for the regional SST measurements not for nearshore or for the SST-differential based upwelling estimate.
</p>
<p>
These last two SST data sets were downloaded from the NOAA ERDDAP server:

<a href="https://coastwatch.pfeg.noaa.gov/erddap/info/ncdcOisst21Agg/index.html">https://coastwatch.pfeg.noaa.gov/erddap/info/ncdcOisst21Agg/index.html</a>

<a href="https://coastwatch.pfeg.noaa.gov/erddap/info/erdAGsstamday/index.html">https://coastwatch.pfeg.noaa.gov/erddap/info/erdAGsstamday/index.html</a>

<a href="https://coastwatch.pfeg.noaa.gov/erddap/info/erdPH2sstamday/index.html">https://coastwatch.pfeg.noaa.gov/erddap/info/erdPH2sstamday/index.html</a>.

<a href="https://coastwatch.pfeg.noaa.gov/erddap/info/esrlIcoads1ge/index.html">https://coastwatch.pfeg.noaa.gov/erddap/info/esrlIcoads1ge/index.html</a> 
</p>
<p>
The SST values were averaged across the thirteen 1
degree by 1 degree boxes which parallel the bathymetry (Figure S1).
</p>
<h3>
References
</h3>
<p>
The AVHRR data were provided by GHRSST and the US National Oceanographic
Data Center. The data were
downloaded from NOAA CoastWatch-West Coast Regional Node and Southwest
Fisheries Science Center’s Environmental Research Division. 
</p>
<p>
The ICOADS data were provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, from their web site at <a href="http://www.esrl.noaa.gov/psd/">http://www.esrl.noaa.gov/psd/</a>
</p>
<p>
Casey KS, Brandon TB, Cornillon P, Evans R (2010). “The past, present,
and future of the AVHRR Pathfinder SST program.” In <em>Oceanography
from space</em>, 273–287. Springer.
</p>
<p>
Simons RA. 2019. ERDDAP. https://coastwatch.pfeg.noaa.gov/erddap . Monterey, CA: NOAA/NMFS/SWFSC/ERD.
</p>
<p>
Walton C, Pichel W, Sapper J, May D (1998). “The development and
operational application of nonlinear algorithms for the measurement of
sea surface temperatures with the NOAA polar-orbiting environmental
satellites.” <em>Journal of Geophysical Research: Oceans</em>,
<b>103</b>(C12), 27999–28012.
</p>
<p>
Huang, Boyin; Liu, Chunying; Banzon, Viva F.; Freeman, Eric; Graham, Garrett; Hankins, Bill; Smith, Thomas M.; Zhang, Huai-Min. (2020): NOAA 0.25-degree Daily Optimum Interpolation Sea Surface Temperature (OISST), Version 2.1. [indicate subset used]. NOAA National Centers for Environmental Information. https://doi.org/10.25921/RE9P-PT57. Accessed November 10, 2020.
</p>

Upwelling Data
==============

<p>
Four upwelling indices were used: a SST nearshore
offshore differential (degree Celcius), the nearshore SST (degree Celcius),  the Ekman Mass Transport
perpendicular to the coast (kg m$^{-1}$ s$^{-1}$), and Ekman Pumping at the tip of India (m s$^{-1}$). The
The SST data were downloaded from the NOAA CoastWatch ERDDAP
server. See the SST data description above. The Ekman Mass Transport and Ekman Pumping were computed from the meridonal and zonal winds from the Reanalysis Data ERA5 monthly Wind velocities downloaded from the Asia Pacific Data-Research Center ERRDAP server. ERA5 uses both the QSAT and ASCAT scatterometer measurements. See Supplement 3 for comparisons of the wind products.
</p>
<p>
The SST differential upwelling indices (degree Celcius) were computed from the ioSST data set (see above in SST section) which is based on AVHRR and thus accurate
close to the coast.The SST-based UPW index is the difference between the neashore box (1 to 5) and a box 3 degrees longitude offshore at the same latitude.
</p>
<p>
The Ekman Mass Transport (EMT) index of upwelling is based upon Ekman’s theory
of mass transport due to wind stress. The index is computed from the 
`ektrx` and `ektry`, which are the x- and y- components of Ekman Transport (kg m$^{-1}$ s$^{-1}$) computed from the longitudinal and latitudinal wind speeds. The functions for computing EMT are below; we used a coast angle of 158 degrees for the India west coast near Kochi (74.5E 11.5N). The `ekrtrx` and `ekrtry` were computed for latitude 8 to 13 and up to 2 degrees longitude from the coast. The values were then averaged to give a coastal average. 

We used monthly winds from the ERA5 Reanalysis product from the European Centre for Medium-Range Weather Forecasts (ECMWF). Using monthly winds (as opposed to daily winds) in the EMT calculation does lead to bias, however monthly winds have been used in many recent papers (see citations in the main text) on upwelling in this region, and our paper sought to test indices that have been previously used. We used the winds closest to the sea surface. See Supplement 3 for a full discussion of the EMT (and Ekman Pumping) calculations along with comparisons of different wind products.

The Reanalysis Data ERA5 monthly 3d Wind velocities were downloaded from </br>

<a href="http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_66d3_10d8_0f3c.html">http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_66d3_10d8_0f3c.html</a>.
 </br> </br>
Our function to compute the EMT perpendicular to the coast was:
</p>
```
getEMT <- function(u, v){
  dat <- list()
  pa <- 1.22 # kg/m3 air pressure
  omega <- 7.272205e-05 #(rad/s)
  f <- 2*omega*sin(pi*dat$latitude/180) #Coriolis parameter
  # Compute tau; get taux and tauy from u and v
  uv_mag <- sqrt(dat$u^2+dat$v^2) # wind speed
  tau  <- pa*Cd(uv_mag)*uv_mag*uv_mag # wind stress
  tauy <- pa*Cd(uv_mag)*uv_mag*dat$v
  taux <- pa*Cd(uv_mag)*uv_mag*dat$u
  EMT <- tau/f
  # MASS TRANSPORT IS PERPENDICULAR TO WIND
  # u/taux positive is west to east wind so into india w coast
  # wind blowing west to east means EMT toward equator (negative)
  EMTy <- -1*taux/f
  # v/tauy negative is wind blowing south toward equator
  # Negative v (south) means EMT to west (off-shore)
  # EMTx is the one that drives upwelling since coast is mostly n-s
  # Negative EMTx = directed offshore = positive upwelling
  EMTx <- tauy/f
  dat$uv_mag <- uv_mag
  dat$tau <- tau
  dat$taux <- taux
  dat$tauy <- tauy
  dat$EMT <- EMT
  dat$ektrx <- EMTx
  dat$ektry <- EMTy
  return(dat)
}
```

```
EMTperp <- function(ektrx, ektry, coast_angle) {
  pi <- 3.1415927
  degtorad <- pi/180.
  alpha <- (180 - coast_angle) * degtorad
  s1 <- cos(alpha)
  t1 <- sin(alpha)
  s2 <- -1 * t1
  t2 <- s1
  perp <- (s1 * ektrx) + (t1 * ektry)
  para <- (s2 * ektrx) + (t2 * ektry)
  return(list(perp=perp, para=para))
}
```

```
getEkmanPump <- function(dat, res){
  # dat is the ektrx and ektry for a grid
  # dat has latitude, longitude, ektrx and ektry
  # res is the grid resolution; 0.25 deg for the ERA5 wind data
  dat$Wey <- NA; dat$Wex <- NA; dat$We <- NA
  lons <- sort(unique(dat$longitude))
  lats <- sort(unique(dat$latitude))

  # constants
  # h is 1 degree lat or lon in meters
  h.meter <- function(lat){
    m_per_deg_lat <- 111132.954 - 559.822 * cos( 2 * pi * lat/180 ) + 
                      1.175 * cos( 4 * pi * lat / 180)
   m_per_deg_lon <- 111132.954 * cos ( pi*lat/180 )
   return(list(lat=m_per_deg_lat, lon=m_per_deg_lon))
  }
  psw <- 1023.6 #kg/m3 density of sea water

  if(length(lats)>2)
    for(lat in (min(lats)+res):(max(lats)-res)){
      dEMTy.lat <- (dat$ektry[dat$latitude==(lat+res)] -
         dat$ektry[dat$latitude==(lat-res)])/
         (2*res*h.meter(lat)$lat)
      Wey <- dEMTy.lat/psw
      dat$Wey[dat$latitude==lat] <- Wey
    }
  if(length(lons)>2)
    for(lon in (min(lons)+res):(max(lons)-res)){
      dEMTx.lon <- (dat$ektrx[dat$longitude==(lon+res)] -
      dat$ektrx[dat$longitude==(lon-res)])/
         (2*res*h.meter(dat$latitude[dat$longitude==(lon+res)])$lon)
      Wex <- dEMTx.lon/psw
      dat$Wex[dat$longitude==lon] <- Wex
    }
  
  dat$We <- dat$Wex+dat$Wey

  return(dat)
}
```

<h3>
References
</h3>
<p>
SST data: These data were provided by GHRSST and the US National
Oceanographic Data Center. This project was supported in part by a grant
from the NOAA Climate Data Record (CDR) Program for satellites. The data
were downloaded from NOAA CoastWatch-West Coast Regional Node and
Southwest Fisheries Science Center’s Environmental Research Division.
</p>
<p>
ERA5: Hersbach, H., Bell, B., Berrisford, P., et al. The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society. 2020; 146: 1999–2049. https://doi.org/10.1002/qj.3803
</p>
<p>
ERA5 website: https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5
<p>
Bakun A (1973). “Coastal upwelling indices, west coast of North
America.” US Department of Commerce, NOAA Technical Report NMFS
SSRF-671.
</p>
<p>
Chamberlain S (2019). rerddap: General Purpose Client for 'ERDDAP' Servers. R package version 0.6.5, https://CRAN.R-project.org/package=rerddap.
</p>
<p>
Mendelssohn R (2020). rerddapXtracto: Extracts Environmental Data from 'ERDDAP' Web Services. R package version 0.4.7, https://CRAN.R-project.org/package=rerddapXtracto.
</p>


Precipitation Data
==================

<p>
We used three precipitation data sets off the southwest coast of India. Two are satellite
derived and one is based on land gauges.
</p>
<p>
The National Climatic Data Center provides basic information on the 
Global Precipitation Climatology Project (GPCP)
Precipitation data set. The data set consists of monthly precipitation
estimates (average mm/day) from January 1979 to the present. The
precipitation estimates merge several satellite and in situ sources into
a final product. Data are provided on a 2.5 degree grid. The GPCP
Precipitation data are provided by the NOAA/NCEI Global Precipitation
Climatology Project and were downloaded from </br>
<a href="https://www.ncei.noaa.gov/data/global-precipitation-climatology-project-gpcp-monthly">https://www.ncei.noaa.gov/data/global-precipitation-climatology-project-gpcp-monthly</a>.
Two boxes were defined, one off the Kerala coast and one off the
Karnataka coast, and the average values of all grid points within these
boxes were used. The boxes are 
Kerala Lat(8.75,  11.25), Lon(73.25,  75.75)
Karnataka Lat(13.75,  16.25), Lon(71.25,  73.75)
</p>
<p>
The land gauge data set is a monthly rainfall (in mm)
area weighted average for each state in India starting from 1901
onwards based on rain gauges. The data are provided by
the India Meteorological Department (Ministry of Earth Sciences). The
1901 to 2014 data were downloaded from the Open Government Data Platform
India <a href="https://data.gov.in">https://data.gov.in</a>. The 2015 and 2016 data were extracted from the yearly
Rainfall Statistics reports (see references).
</p>
<p>
NASA’s Tropical Rainfall Measuring Mission (TRMM) website provides
background on the TRMM precipitation data
(<a href="https://pmm.nasa.gov/">https://pmm.nasa.gov/</a>). 1997 to
2015 monthly precipitation estimates on a 0.25 degree grid were
downloaded from the Tropical Rainfall Measuring Mission (TRMM) website.
The data were averaged in the 2.5 x 2.5 degree boxes 1 to 13 used for the other
satellite data.
</p>
<h3>
References
</h3>
<p>
Adler R, Huffman G, Chang A, Ferraro R, Xie P, Janowiak J, Rudolf B,
Schneider U, Curtis S, Bolvin D, Gruber A, Susskind J, Arkin P (2003).
“The Version 2 Global Precipitation Climatology Project (GPCP) Monthly
Precipitation Analysis (1979-Present).” <em>Journal of
Hydrometeorology</em>, <b>4</b>, 1147-1167.
</p>
<p>
Adler R, Wang J, Sapiano M, Huffman G, Chiu L, Xie PP, Ferraro R,
Schneider U, Becker A, Bolvin D, Nelkin E, Gu G, Program NC (2016).
“Global Precipitation Climatology Project (GPCP) Climate Data Record
(CDR), Version 2.3 (Monthly).” National Centers for Environmental
Information. doi:
<a href="https://doi.org/10.7289/V56971M6">10.7289/V56971M6</a>.
</p>
<p>
Purohit MK, Kaur S (2016). “Rainfall Statistics of India - 2016.” India
Meterorological Department (Ministry of Earth Sciences).
<a href="http://hydro.imd.gov.in/hydrometweb/">http://hydro.imd.gov.in/hydrometweb/</a>.
</p>
<p>
Kothawale DR, Rajeevan M (2017). Monthly, seasonal and annual rainfall time series for all-India, homogeneous regions and meteorological subdivisions: 1871-2016. Report, Indian Institute of Tropical Meteorology..
</p>
<p>
NCEI (2017). Global Precipitation Climatology Project Monthly Product Version 2.3. Retrieved from National Centers for Environmental Information website:  
<a href="https://www.ncei.noaa.gov/data/global-precipitation-climatology-project-gpcp-monthly/access/">https://www.ncei.noaa.gov/data/global-precipitation-climatology-project-gpcp-monthly/access/</a>.
</p>

<hr />

Chlorophyll Data
================

<p>
We used Chlorophyll-a products developed by the Ocean Biology Processing
Group in the Ocean Ecology Laboratory at the NASA Goddard Space Flight
Center.
</p>
<p>
For 1997 to 2002, we used the Chlorophyll-a 2014.0 Reprocessing
(R2014.0) product from the Sea-viewing Wide Field-of-view Sensor
(SeaWiFS) on the Orbview-2 satellite. These data are on a 0.1 degree
grid with units of mg m$^{-3}$. See reference below.
</p>
<p>
For 2003 to 2017, we used the MODIS-Aqua product on a 4km grid  with units of mg m$^{-3}$. These
CHL data are taken from measurements gathered by the Moderate Resolution
Imaging Spectroradiometer (MODIS) on NASA’s Aqua Spacecraft. See
reference below.
</p>
<p>
Both CHL data sets were downloaded from the NOAA ERDDAP server:

<a href="https://coastwatch.pfeg.noaa.gov/erddap/info/erdSW1chlamday/index.html">https://coastwatch.pfeg.noaa.gov/erddap/info/erdSW1chlamday/index.html</a>

<a href="https://coastwatch.pfeg.noaa.gov/erddap/info/erdMH1chlamday/index.html">https://coastwatch.pfeg.noaa.gov/erddap/info/erdMH1chlamday/index.html</a>.
</p>
<p>
The CHL values were averaged across the thirteen 1
degree by 1 degree boxes which parallel the bathymetry (Figure S1).
</p>
<h3>
References
</h3>
<p>
NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean
Biology Processing Group; (2014): SeaWiFS Ocean Color Data; NASA Goddard
Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing
Group.
<a href="https://dx.doi.org/10.5067/ORBVIEW-2/SEAWIFS_OC.2014.0">https://dx.doi.org/10.5067/ORBVIEW-2/SEAWIFS\_OC.2014.0</a>
</p>
<p>
NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean
Biology Processing Group. Moderate-resolution Imaging Spectroradiometer
(MODIS) Aqua Chlorophyll Data; 2014 Reprocessing. NASA OB.DAAC,
Greenbelt, MD, USA.
<a href="https://dx.doi.org/10.5067/AQUA/MODIS/L3M/CHL/2014">https://dx.doi.org/10.5067/AQUA/MODIS/L3M/CHL/2014</a>
</p>
<p>
Hu C, Lee Z, Franz B (2012). “Chlorophyll aalgorithms for oligotrophic
oceans: A novel approach based on three-band reflectance difference.”
<em>Journal of Geophysical Research: Oceans</em>, <b>117</b>(C1).
</p>

Ocean Climate Indices
=========

<p>
We used the following ocean climate indices: Oceanic Nino Index, Dipole Mode Index, Pacific Decadal Oscillation index and Atlantic Multidecadal Oscillation index.
</p>
<p>
The Oceanic Nino Index (ONI) is one of the primary indices used to monitor the El Nino-Southern Oscillation (ENSO). The ONI index is 3 month running mean of ERSST.v5 SST anomalies in the
Niño 3.4 region (5°N-5°S, 120°-170°W)\], based on centered 30-year base
periods updated every 5 years.  
The ONI was downloaded from the National Weather Service Climate Prediction Center 

<a href="http://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt">http://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt</a>
<p>
The DMI is the monthly Dipole Mode Index. The DMI is defined by the SSTA (SST anomaly) difference between the
western Indian Ocean (10°S–10°N, 50°E–70°E) and the southeastern Indian
Ocean (10°S–0°, 90°E–110°E). The data were downloaded from

<a href="https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/dmi.long.data">https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/dmi.long.data</a> The original data source is the
Japan Agency for Marine-Earth Science and Technology (JAMSTEC) via the page 

<a href="http://www.jamstec.go.jp/frcgc/research/d1/iod/e/iod/about_iod.html">http://www.jamstec.go.jp/frcgc/research/d1/iod/e/iod/about_iod.html</a>.
</p>
<p>
The PDO is the Pacific Decadal Oscillation. The PDO index is defined by the sea surface temperature anomaly over the North Pacific Ocean (north of 20°N). The data were downloaded from the NOAA Physical Sciences Laboratory:

<a href="https://psl.noaa.gov/tmp/gcos_wgsp/data.143.131.2.6.325.11.4.55">https://psl.noaa.gov/tmp/gcos_wgsp/data.143.131.2.6.325.11.4.55</a>
</p>
<p>
The AMO is the Atlantic Multidecadal Oscillation. The AMO index is based on sea-surface temperature anomalies in the North Atlantic Ocean. The data were downloaded from the NOAA Physical Sciences Laboratory:

<a href="https://psl.noaa.gov/tmp/gcos_wgsp/data.143.131.2.6.325.11.25.33">https://psl.noaa.gov/tmp/gcos_wgsp/data.143.131.2.6.325.11.25.33</a>
</p>
<h3>
References
</h3>
<p>
Saji NH, Yamagata T (2003). “Possible impacts of Indian Ocean Dipole
mode events on global climate.” <em>Climate Research</em>, <b>25</b>(2),
151-169.
</p>
<p>Climate Prediction Center. El Nino-Southern Oscillation. 

<a href="https://www.cpc.ncep.noaa.gov/products/precip/CWlink/MJO/enso.shtml">
https://www.cpc.ncep.noaa.gov/products/precip/CWlink/MJO/enso.shtml</a>
</p>
<p>NOAA Physical Sciences Laboratory. Pacific Decadal Oscillation. 

<a href="https://psl.noaa.gov/pdo/">
https://psl.noaa.gov/pdo/</a>
</p>
<p>
van den Dool HM, Saha S, Johansson Å (2000). “Empirical orthogonal teleconnections.” <em>Journal of Climate</em>, <b>13</b>, 1421–1435.
</p>
<p>
Mantua NJ, Hare SR, Zhang Y, Wallace JM, Francis RC (1997). “A Pacific interdecadal climate oscillation with impacts on salmon production.” <em>Bulletin of the American Meteorology Society</em>, <b>78</b>, 1069-1079.
</p>
<p>
Enfield DB, Mestas-NuÃ±ez AM, Trimble PJ (2001). “The Atlantic Multidecadal Oscillation and its relation to rainfall and river flows in the continental U.S.” <em>Geophysical Research Letters</em>, <b>28</b>, 2077–2080. 
</p>
<hr />

