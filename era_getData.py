import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type':'reanalysis',
        'format':'netcdf',
	"area": "47/9/46/10",
        'variable':[
            'geopotential','temperature','u_component_of_wind',
            'v_component_of_wind', 'relative_humidity'
        ],
        'pressure_level':[
            '400','500','600','625','650','675','700','725', '750', '775','800','825','850',
            '875', '1000'
        ],
        'year':'2015',
        'month':[
            '01'
        ],
        'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31'
        ],
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ],
        'format':'netcdf'
    },
    'plevel.nc')


import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'variable':['geopotential', '2m_dewpoint_temperature', 'surface_thermal_radiation_downwards', 'surface_solar_radiation_downwards',
        'Total precipitation','2m_temperature', 'TOA incident solar radiation',
            'friction_velocity','instantaneous_moisture_flux','instantaneous_surface_sensible_heat_flux'
        ],
        'product_type':'reanalysis',
		"area": "47/9/46/10",
        'year':'2015',
        'month':[
            '01'
        ],
        'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31'
        ],
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ],
        'format':'netcdf'

    },
    'surface.nc')
