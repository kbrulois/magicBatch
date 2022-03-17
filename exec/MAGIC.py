#!/usr/local/bin/python3

import os
import sys
import argparse
import numpy as np 
import datatable as dt
import magicBatch

class ArgumentParserError(Exception):
	pass

class NewArgumentParser(argparse.ArgumentParser):
	def error(self, message):
		print(message)
		sys.exit(0)


def parse_args(args):
	p = NewArgumentParser(description='run MAGIC')
	a = p.add_argument_group('data loading parameters')
	a.add_argument('-d', '--data', metavar='D',
					help='File path of input data csv file (rows represent cells, columns represent features).')
	a.add_argument('-m', '--mar_mat_in', metavar='M', 
					help='File path of input data csv file (rows represent cells, columns represent features) to be used directly in Marcov matrix calculation.')
	a.add_argument('-o1', '--mar_mat_out', metavar='O1',
				   help='File path of where to save the Marcov matrix in sparse coo format (in csv format).')
	a.add_argument('-o2', '--dif_map_out', metavar='O2',
				   help='File path of where to save the diffusion map data (in csv format).')
	a.add_argument('-o3', '--data_out', metavar='O3',
				   help='File path of where to save the imputed data (in csv format).')

	m = p.add_argument_group('MAGIC parameters')
	m.add_argument('-t', metavar='T',
				   help='t paramter to use when running MAGIC.')
	m.add_argument('-p', '--pca-components', metavar='P', default=20, type=int,
				   help='Number of pca components to use when running MAGIC (Default = 20).')
	m.add_argument('--pca-non-random', default=True, action='store_false',
				    help='Do not used randomized solver in PCA computation.')
	m.add_argument('-k', metavar='K', default=9, type=int,
					help='Number of nearest neighbors to use when running MAGIC (Default = 9).')
	m.add_argument('-ka', metavar='KA', default=3, type=int,
					help='knn-autotune parameter for running MAGIC (Default = 3).')
	m.add_argument('-e', '--epsilon', metavar='E', default=1, type=int,
					help='Epsilon parameter for running MAGIC (Default = 1).')
	m.add_argument('-r', '--rescale', metavar='R', type=float, default=0,
					help='Percentile to rescale data to after running MAGIC (Default = 99).')
	m.add_argument('-rm', '--rescale_method', metavar='RM', default='classic',
					help='Method to use for rescaling data after running MAGIC (Default = classic).')

	w = p.add_argument_group('Diffusion Map parameters')
	w.add_argument('-c', '--n_diffusion_components', metavar='C', default=10, type=int,
					help='Number of diffusion map components to calculate (Default = 10).')

	try:
		return p.parse_args(args)
	except ArgumentParserError:
		raise


def main(args: list = None):
	args = parse_args(args)
	print(args)
	try:
		magicBatch.MAGIC_core.magicBatch(data=args.data, mar_mat_input=args.mar_mat_in, n_pca_components=args.pca_components, random_pca=args.pca_non_random, t=args.t, k=args.k, ka=args.ka, epsilon=args.epsilon, n_diffusion_components=args.n_diffusion_components, rescale_percent = args.rescale, rescale_method = args.rescale_method, csv_l=args.mar_mat_out, csv_d=args.dif_map_out, csv_i=args.data_out)

	except:
		raise

if __name__ == '__main__':
	main(sys.argv[1:])
