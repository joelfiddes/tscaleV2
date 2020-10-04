# tscale_lib

"""
python2

Core tscale methods

Example:


Vars:


Details:


"""


def tscale_1d(tscaleSrc, wdir, startDate, endDate, dataset, windCor, plapse):
	'''
	Make listpoints

	Runs: makeListpoints2_points

	Args:
		tmappSrc:
		wdir:
	'''

	cmd = [
			"python",
			tscaleSrc + "/toposcale/tscale_run.py",
			wdir,
			home,
			startDate,
			endDate,
			dataset,
			windCor,
			plapse
		]
	subprocess.check_output(cmd)

def tscale_3d(tscaleSrc, wdir, runmode, startDate, endDate):
	'''
	Make listpoints

	Runs: makeListpoints2_points

	Args:
		tmappSrc:
		wdir:
	'''
	cmd = [
			"python",
			tscaleSrc + "/toposcale/tscale3D.py",
			wdir,
			runmode,
			startDate,
			endDate,
			"HRES",
			'1'

		]

	subprocess.check_output(cmd)
