# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 11:17:40 2022

@author: t
"""


import xarray as xr
from scipy.optimize import curve_fit
from shapely.geometry.polygon import Polygon
import shapely.ops as ops
from functools import partial
from pyproj import Proj, transform
import numpy as np
from pathlib import Path


def cal_area(lat, lon, res=0.008333 / 2):
	lat_m = lat
	lon_m = lon
	if lat_m > 89:
		lat_m = lat_m - res
	if lat_m < -89:
		lat_m = lat_m + res
	geom = Polygon(
		[(lon_m + res, lat_m + res),
		 (lon_m + res, lat_m - res),
		 (lon_m - res, lat_m - res),
		 (lon_m - res, lat_m + res)])
	# print(geom)
	geom_area = ops.transform(
		partial(
			transform,
			Proj(init='EPSG:4326'),
			Proj(
				proj='aea',
				lat_1=geom.bounds[1],
				lat_2=geom.bounds[3]
			)
		),
		geom)
	grid_area = geom_area.area
	# print('grid_area',grid_area,'m^2')
	return grid_area


def get_t2(da_t2, lat, lon):
	# if lon<0:
	#     lon_360=lon+360
	# else:
	#     lon_360=lon
	da_t2 = da_t2.interp(longitude=lon, latitude=lat)
	'''change lontitude from -180-180 to 0-360'''

	t2 = da_t2.t2m.data - 273.5
	return t2


def cal_area(lat, lon, res=0.008333 / 2):
	lat_m = lat
	lon_m = lon
	if lat_m > 89:
		lat_m = lat_m - res
	if lat_m < -89:
		lat_m = lat_m + res
	geom = Polygon(
		[(lon_m + res, lat_m + res),
		 (lon_m + res, lat_m - res),
		 (lon_m - res, lat_m - res),
		 (lon_m - res, lat_m + res)])
	# print(geom)
	geom_area = ops.transform(
		partial(

			transform,
			Proj(init='EPSG:4326'),
			Proj(
				proj='aea',
				lat_1=geom.bounds[1],
				lat_2=geom.bounds[3]
			)
		),
		geom)
	grid_area = geom_area.area
	# print('grid_area',grid_area,'m^2')
	return grid_area  # 'm^2'


def get_pop_d(da_pop, lat, lon):
	pop = da_pop.interp(x=lon, y=lat).data
	area = cal_area(lat, lon, )
	pop_d = (pop) / area * 10000

	return pop_d  # poplation density （cap/ha-1）


def get_ahf(ds_ahf, lat, lon):
	ahf = ds_ahf.interp(x=lon, y=lat).to_array()
	return ahf / 10 ** 5


def fit_ahf(t2m, ahf, pop):
	def cal_ahf(t, a_f0, a_f1, a_f2):
		hdd_t = 16  # cdd threshold
		cdd_t = 26  # hdd threshold
		cdd = np.where(t - cdd_t > 0, t - cdd_t, 0)

		hdd = np.where(hdd_t - t > 0, hdd_t - t, 0)

		return pop * (a_f0 + a_f1 * cdd + a_f2 * hdd)

	coefficient, time = curve_fit(cal_ahf, t2m, ahf)

	ahf_fitted = cal_ahf(t2m, coefficient[0], coefficient[1], coefficient[2])
	# print('af_base,af_cdd,af_hdd',coefficient)
	return coefficient, ahf_fitted


def cal_coef(da_t2, da_ahf, da_pop, lat, lon):
	# read data
	t2 = get_t2(da_t2, lat, lon)  # 12 monthly mean temperature
	# ahf=get_ahf(ds_ahf,lat,lon) # 12 monthly mean anthropologic heat flux
	ahf = da_ahf.data[:, 0] / 10 ** 5
	pop_d = get_pop_d(da_pop, lat, lon)  # poplation density （cap/ha-1）
	# curve_fit
	c, f = fit_ahf(t2, ahf, pop_d)
	return c, t2, pop_d, ahf


def derive_hours_ahf(coef, t2, pop_d, df_weight, hours, tz=0):
	hdd_t = 16
	cdd_t = 26
	cdd = np.where(t2 - cdd_t > 0, t2 - cdd_t, 0)
	hdd = np.where(hdd_t - t2 > 0, hdd_t - t2, 0)
	ahf_day = pop_d * (coef[0] + coef[1] * cdd + coef[2] * hdd)
	''''
	tz:timezone
	hours: UTC hours
	'''
	hours = (hours + tz) % 24
	if t2 < 12.4:
		w = df_weight.iloc[hours, 0]
		if t2 < 16.95:
			w = df_weight.iloc[hours, 1]
			if t2 < 20.95:
				w = df_weight.iloc[hours, 2]
			else:
				w = df_weight.iloc[hours, 2]
	ahf_hour = ahf_day * w
	return ahf_hour


if __name__ == '__main__':
	import pandas as pd

	# df_city=pd.read_excel('../data/citylist.xls')
	lat = 71.3039629414
	lon = -156.730205754

	''' temperature from era5 
		../data/download_era5.py
	'''
	fn_t2 = '../data/era5/reanalysis-era5-single-levels-monthly-means.nc'
	da_t2 = xr.open_dataset(fn_t2)
	da_t2 = da_t2.assign_coords(
		{"longitude": (xr.where(da_t2.longitude > 180, da_t2.longitude - 360, da_t2.longitude))})
	# da_t2=da_t2.rename_dims({"longitude": 'lon','latitude':'lat'})

	''' 
	population  from LandScan Global 2013 

	for more detail ../data/LandScan Global 2013/LandScan 2013 Metadata.htm
	'''
	fn_pop = "../data/LandScan Global 2013/lspop2013/w001001.adf"
	da_pop = xr.open_rasterio(fn_pop)

	''' 
	ahf  weight for derive hourly data
	'''
	fn_weight = '../data/weight.csv'
	df_weight = pd.read_csv(fn_weight, ).iloc[:, 1:]

	# 12 monthly mean anthropologic heat flux
	''' 
	AHF (Q_F) from AH4GUC ([Varquez, 2021]
	(https://urbanclimate.tse.ens.titech.ac.jp/database/AHE/AH4GUC/present/))
	'''
	ls_fn_ahf = list(Path('../data/ahf').glob("*mean_UTC.tif"))[:-1]
	da_ahf = xr.concat([xr.open_rasterio(fn_ahf).interp(x=lon, y=lat) for fn_ahf in ls_fn_ahf], dim='mon')
	# fn_ahf='../data/ahf_m.nc'
	# ds_ahf=xr.open_dataset(fn_ahf)
	'''get coefficient'''
	coef, t2, pop_d, ahf = cal_coef(da_t2, da_ahf, da_pop, lat, lon)

	ahf_hour = derive_hours_ahf(coef, pop_d=pop_d, df_weight=df_weight, hours=15, t2=t2[5])
	print(ahf_hour)