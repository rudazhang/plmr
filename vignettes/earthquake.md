# Earthquake Data

The Advanced National Seismic System (ANSS) is a collaboration of the U.S. Geological Survey (USGS)
and regional, state, and academic partners that collect and analyze data on significant earthquakes
to provide near real-time (generally within 10 to 30 minutes) information.

The ANSS Comprehensive Earthquake Catalog (ComCat) contains earthquake source parameters,
e.g. hypocenters, magnitudes, phase picks and amplitudes,
produced by contributing seismic/seismograph networks.

[ComCat Documentation](https://earthquake.usgs.gov/data/comcat/):
[catalogs](https://earthquake.usgs.gov/data/comcat/catalog/);
[contributors](https://earthquake.usgs.gov/fdsnws/event/1/contributors);

[Data Availability](https://earthquake.usgs.gov/data/comcat/data-availability.php)
Earthquakes occurring outside the US and smaller than about magnitude 4.5
can be difficult for the USGS to locate if there are not enough data.
The USGS continues to receive data from observatories throughout the world
for several months after the events occur.
Using those data, we add new events and revise existing events in later publications.

In most cases, we locate and report on earthquakes worldwide M5.0+ within 30 min.

[Seismicity Map 1900-2013](https://www.usgs.gov/natural-hazards/earthquake-hazards/information-region)


## API

[Earthquake Catalog API Documentation](https://earthquake.usgs.gov/fdsnws/event/1/)

URL pattern: `<root>/<method>?<parameters>`

URL root: https://earthquake.usgs.gov/fdsnws/event/1/

| API method | Service                 |
|------------|-------------------------|
| `count`    | Count number of returns |
| `query`    | Submit a data request   |

Parameters:
- `format`: output format; "csv", "geojson";
- `limit`: limit of the number of returned events, defaults to 20000 and cannot be increased;
- location parameters: default to full range;
- `starttime`, `endtime`: range of event time; "2019-11-19";
- `mindepth`, `maxdepth`: event depth range in km, defaults to [-100, 1000] and cannot be extended;
- `minmagnitude`, `maxmagnitude`: event magnitude range, no default;
- `eventtype`: type of seismic event; "earthquake";

## Data Dictionary

`net`: identifier of preferred contributor;
`id`: event identifiers specific to a data center;
`status`: review status; "reviewed", "automatic", "deleted";
automatic events are directly posted by automatic processing systems,
reviewed events have been looked at by a human;
`types`: a comma-separated list of product types associated to this event;
moment-tensor, focal-mechanism, shakemap, losspager, dyfi;

`type`: earthquake, quarry, volcanic eruption, nuclear explosion, etc.

Coordinates of the epicenter are given in the WGS84 reference frame.
Reference elevation depends on the method a seismic network used to locate the earthquake,
may be the WGS84 geoid, mean sea-level, or the average elevation of the seismic stations
which provided arrival-time data for the earthquake location.
- `latitude`: [-90, 90], degree;
- `longitude`: [-180, 180], degree;
- `depth`: [0, 1000], kilometer;
  defaults to 10 km for mid-continental areas & mid-ocean ridges, 33/35 km for shallow earthquakes;
- `place`: description of named geographic region near the event;

Location errors:
The principal errors are the major axes of the error ellipsoid.
"Unknown" if the contributing seismic network does not supply uncertainty estimates.
- `horizontalError`: horizontal location error, kilometer;
  length of the largest projection of the three principal errors on a horizontal plane.
- `depthError`: vertical location error, kilometer;
  the largest projection of the three principal errors on a vertical line.
- `gap`: largest azimuthal gap between azimuthally adjacent stations, degrees;
  large azimuthal gap implies large location and depth uncertainties;
- `dmin`: horizontal distance from the epicenter to the nearest station, degree;
  large horizontal distance implies large depth uncertainty;
- `rms`: root-mean-square travel time residual, sec;
  smaller numbers reflect a better fit of the data to determin location.

All times use ISO8601 Date/Time format, defaults to UTC.
- `time`: time when the event initiates rupture;
- `updated`: time when the event was most recently updated;

Magnitude:
`mag`: USGS official magnitude for the event;
`magError`: estimated standard error of the magnitude;
`magType`: method used to calculate the magnitude;
The various types of magnitude measure different aspects of the seismic radiation,
e.g., low-frequency energy vs. high-frequency energy.
[Magnitude types](https://www.usgs.gov/natural-hazards/earthquake-hazards/science/magnitude-types):
"Mw"/"Mww", W-phase; "Mwc", long-period surface waves; "Mwb", long-period body-waves;
"Mwr", whole seismogram at regional distances; "mb", short-period body wave; etc.

`locationSource`: network that originally authored the reported location;
`nst`: number of seismic stations reported P- and S-arrival times to determine location of the event;
`magSource`: network that originally authored the reported magnitude;
`magNst`: number of seismic stations used to calculate the magnitude;

