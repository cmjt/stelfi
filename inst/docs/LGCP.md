    library(stelfi)

NZ murders
==========

The Data
--------

    data(murders_nz)
    dim(murders_nz)

\[1\] 967 13

    head(murders_nz)

Latitude Longitude Sex Age Date Year Cause Killer 1 -43.63394 171.6442
Male 41 Jan 5 2004 stabbing friend 2 -43.28563 172.1305 Male 46 Jan 8
2004 pick axe wounds friend 3 -36.92575 174.8498 Male 0 Jan 15 2004
asphyxiation (suffocation) mother 4 -43.55006 172.6327 Female 46 Feb 1
2004 blunt force trauma partner 5 -40.73297 175.1195 Male 10 Feb 2 2004
stabbing father 6 -40.73273 175.1193 Female 2 Feb 2 2004 stabbing father
Name Full\_date Month Cause\_cat Region 1 Donald Linwood 2004-01-05
January Violent weapon Canterbury 2 James Weeks 2004-01-08 January
Violent weapon Canterbury 3 Gabriel Harrison-Taylor 2004-01-15 January
Asphyxia Auckland 4 Odette Lloyd-Rangiuia 2004-02-01 February Blunt
force trauma Canterbury 5 Te Hau OCarroll 2004-02-02 February Violent
weapon Wellington 6 Ngamata OCarroll 2004-02-02 February Violent weapon
Wellington

    ## Warning in kable_styling(., bootstrap_options = "striped", full_width = FALSE):
    ## Please specify format in kable. kableExtra can customize either HTML or LaTeX
    ## outputs. See https://haozhu233.github.io/kableExtra/ for details.

<table>
<caption>Number of murders by category</caption>
<thead>
<tr class="header">
<th align="left">Cause</th>
<th align="right">Number</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Asphyxia</td>
<td align="right">47</td>
</tr>
<tr class="even">
<td align="left">Blunt force trauma</td>
<td align="right">285</td>
</tr>
<tr class="odd">
<td align="left">Car crash</td>
<td align="right">112</td>
</tr>
<tr class="even">
<td align="left">Drowning</td>
<td align="right">16</td>
</tr>
<tr class="odd">
<td align="left">Drugs</td>
<td align="right">10</td>
</tr>
<tr class="even">
<td align="left">Fire</td>
<td align="right">16</td>
</tr>
<tr class="odd">
<td align="left">Other</td>
<td align="right">72</td>
</tr>
<tr class="even">
<td align="left">Violent weapon</td>
<td align="right">409</td>
</tr>
</tbody>
</table>

    ## Warning in kable_styling(., bootstrap_options = "striped", full_width = FALSE):
    ## Please specify format in kable. kableExtra can customize either HTML or LaTeX
    ## outputs. See https://haozhu233.github.io/kableExtra/ for details.

<table>
<caption>Number of murders by year</caption>
<thead>
<tr class="header">
<th align="left">Year</th>
<th align="right">Number</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">2004</td>
<td align="right">38</td>
</tr>
<tr class="even">
<td align="left">2005</td>
<td align="right">72</td>
</tr>
<tr class="odd">
<td align="left">2006</td>
<td align="right">66</td>
</tr>
<tr class="even">
<td align="left">2007</td>
<td align="right">55</td>
</tr>
<tr class="odd">
<td align="left">2008</td>
<td align="right">71</td>
</tr>
<tr class="even">
<td align="left">2009</td>
<td align="right">94</td>
</tr>
<tr class="odd">
<td align="left">2010</td>
<td align="right">71</td>
</tr>
<tr class="even">
<td align="left">2011</td>
<td align="right">62</td>
</tr>
<tr class="odd">
<td align="left">2012</td>
<td align="right">64</td>
</tr>
<tr class="even">
<td align="left">2013</td>
<td align="right">59</td>
</tr>
<tr class="odd">
<td align="left">2014</td>
<td align="right">55</td>
</tr>
<tr class="even">
<td align="left">2015</td>
<td align="right">65</td>
</tr>
<tr class="odd">
<td align="left">2016</td>
<td align="right">56</td>
</tr>
<tr class="even">
<td align="left">2017</td>
<td align="right">45</td>
</tr>
<tr class="odd">
<td align="left">2018</td>
<td align="right">76</td>
</tr>
<tr class="even">
<td align="left">2019</td>
<td align="right">18</td>
</tr>
</tbody>
</table>

    ## Warning in kable_styling(., bootstrap_options = "striped", full_width = FALSE):
    ## Please specify format in kable. kableExtra can customize either HTML or LaTeX
    ## outputs. See https://haozhu233.github.io/kableExtra/ for details.

<table>
<caption>Number of murders by province</caption>
<thead>
<tr class="header">
<th align="left">Region</th>
<th align="right">Number</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Auckland</td>
<td align="right">267</td>
</tr>
<tr class="even">
<td align="left">Bay of Plenty</td>
<td align="right">94</td>
</tr>
<tr class="odd">
<td align="left">Canterbury</td>
<td align="right">116</td>
</tr>
<tr class="even">
<td align="left">Gisborne</td>
<td align="right">17</td>
</tr>
<tr class="odd">
<td align="left">Hawke's Bay</td>
<td align="right">44</td>
</tr>
<tr class="even">
<td align="left">Manawatu-Wanganui</td>
<td align="right">65</td>
</tr>
<tr class="odd">
<td align="left">Marlborough</td>
<td align="right">6</td>
</tr>
<tr class="even">
<td align="left">Nelson</td>
<td align="right">16</td>
</tr>
<tr class="odd">
<td align="left">Northland</td>
<td align="right">49</td>
</tr>
<tr class="even">
<td align="left">Otago</td>
<td align="right">30</td>
</tr>
<tr class="odd">
<td align="left">Southland</td>
<td align="right">19</td>
</tr>
<tr class="even">
<td align="left">Taranaki</td>
<td align="right">27</td>
</tr>
<tr class="odd">
<td align="left">Tasman</td>
<td align="right">6</td>
</tr>
<tr class="even">
<td align="left">Waikato</td>
<td align="right">108</td>
</tr>
<tr class="odd">
<td align="left">Wellington</td>
<td align="right">89</td>
</tr>
<tr class="even">
<td align="left">West Coast</td>
<td align="right">11</td>
</tr>
</tbody>
</table>

    data(nz) ## SpatialPolygonsDataFrame of NZ (NZTM projection)
    area_nz <- raster::area(nz) ## (m)
    area_nzkm2 <- area_nz/1000^2 ## according to Google NZ is 268,021 km2
    spatial_murder_rate <- nrow(murders_nz)/area_nzkm2 
    temporal_murder_rate <- nrow(murders_nz)/length(table(murders_nz$Year)) 
    st_murder_rate <- (nrow(murders_nz)/area_nzkm2)/length(table(murders_nz$Year)) 

Rate of murders per km**<sup>2</sup> across NZ is calculated as 0.19,
0.072, 0.022, 0.116, 0.076, 0.044, 0.075, 0.03, 0.029, 0.121, 0.039,
0.12, 0.041; there are roughly 60.438. The spatio-temporal murder rate
across NZ 2004--2019 is 0.012, 0.004, 0.001, 0.007, 0.005, 0.003, 0.005,
0.002, 0.002, 0.008, 0.002, 0.007, 0.003 (rate per km**<sup>2</sup> per
year).

### Transform to NZTM

Transform `data.frame` to `SpatialPointsDataFrame`

    murders_sp <- murders_nz
    ## project longitude & latitude to NZTMs
    coordinates(murders_sp) <- c("Longitude","Latitude")
    proj4string(murders_sp) <- CRS("+proj=longlat +datum=WGS84")
    murders_sp <-  spTransform(murders_sp, 
                               CRS("+proj=nzmg +lat_0=-41.0 +lon_0=173.0 +x_0=2510000.0 +y_0=6023150.0 +ellps=intl +units=m"))

![Locations of recorded (n = 967) murders in NZ
2004--2019](LGCP_files/figure-markdown_strict/plot-1.png)

log-Gaussian Cox process
------------------------

### Using INLA

**Steps below closely follow this [INLA-SPDE
tutorial](https://becarioprecario.bitbucket.io/spde-gitbook/ch-lcox.html)**.

#### Creating the mesh

Typically when analysing point pattern data the point locations are not
specified as the mesh nodes (i.e., locations are not given as an
argument to `inla.mesh.2d()`). Instad we can supply the coordinates of
the point pattern window (domain).

    mesh <- inla.mesh.2d(loc.domain = coordinates(nz) ,
                         max.edge = c(86000, 100000), cutoff = 5000)
