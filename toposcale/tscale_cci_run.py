# https://www.evernote.com/Home.action#n=cf913e83-0c04-469c-ae65-fab96083d751&s=s500&ses=4&sh=5&sds=2&x=addresses&


import tscale_cci
import datetime
import pandas as pd
from tqdm import tqdm
startYear=1980
endYear=2018#2001
coords_file = "/cluster/home/fiddes/nn9606k/tscale/coordinates.dat"
era_dir = "/cluster/home/fiddes/nn9606k/era5"
out_dir	= "/cluster/home/fiddes/nn9606k/tscale"


years= range(startYear,endYear+1)

for year in tqdm(years):
        print(year)
        startPeriods = pd.date_range(str(year)+'-01-01', periods=46, freq='8d')

        for startIndex in tqdm(range(0,startPeriods.size)):

                start=startPeriods[startIndex]
                # always start last period at 23rd Dec to ensure 8 day period till 31 Dec
                # python indexing (starts 0) so startIndex 45 == period 46
                if startIndex==45:
                        start = str(year)+'-12-23'
                        t= pd.to_datetime(start)
                        end = t+datetime.timedelta(days=8)

                if startIndex<45:
                        t= pd.to_datetime(start)
                        end = t+datetime.timedelta(days=9)



                tscale_cci.main( coords_file,era_dir , out_dir, str(start),  str(end), startIndex )



