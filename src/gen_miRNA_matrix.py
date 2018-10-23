# copyright: yueshi@usc.edu
import pandas as pd 
import hashlib
import os
import numpy as np
import numpy
from utils import logger
from numpy import nan
def file_as_bytes(file):
    with file:
        return file.read()

def extractMatrix(dirname):
	'''
	return a dataframe of the miRNA matrix, each row is the miRNA counts for a file_id

	'''
	count = 0

	miRNA_data = []
	for idname in os.listdir(dirname):
		# list all the ids 
		if idname.find("-") != -1:
			idpath = dirname +"/" + idname

			# all the files in each id directory
			for filename in os.listdir(idpath):
				# check the miRNA file
				if filename.find("-") != -1:

					filepath = idpath + "/" + filename
					df = pd.read_csv(filepath,sep="\t")
					# columns = ["miRNA_ID", "read_count"]
					if count ==0:
						# get the miRNA_IDs 
						miRNA_IDs = df.miRNA_ID.values.tolist()

					id_miRNA_read_counts = [idname] + df.read_count.values.tolist()
					miRNA_data.append(id_miRNA_read_counts)


					count +=1
					# print (df)
	columns = ["file_id"] + miRNA_IDs
	df = pd.DataFrame(miRNA_data, columns=columns)
	return df

def extractLabel(inputfile, casefile):
	df = pd.read_csv(inputfile, sep="\t")
	df_case = pd.read_csv(casefile, sep="\t")
	#
	# print (df[columns])
	df['label'] = df['cases.0.samples.0.sample_type']

	df = pd.read_csv("/Users/xin/Downloads/EE542/Labs/lab10/data/files_meta.tsv", sep="\t")
	df_case = pd.read_csv("/Users/xin/Downloads/EE542/Labs/lab10/data/cases_meta.tsv", sep="\t")
	df_map = pd.read_csv("/Users/xin/Downloads/EE542/Labs/lab10/data/file_case_id_DNA.csv")

	df.loc[df['cases.0.samples.0.sample_type'].str.contains("Normal"), 'label'] = 0
	df.loc[~df['cases.0.samples.0.sample_type'].str.contains("Normal"), 'label'] = 1

	for i in range(df.shape[0]):
		cur = df_map.loc[df_map['file_id'] == df.loc[i, 'file_id']]
		case_id = cur['case_id'].item()
		cur_case = df_case.loc[df_case['case_id'] == case_id]
		cur_label = cur_case.primary_site.item()
		if df.loc[i,'label'] != 0:
			df.loc[i,'label'] = cur_label

	labels = {0: 0, 'Colon': 1, 'Breast': 2, 'Esophagus': 3, nan: 60, 'Kidney': 4, 'Prostate gland': 5,
			  'Brain': 6, 'Adrenal gland': 7, 'Thymus': 8, 'Lymph nodes': 9,
			  'Base of tongue': 10, 'Other and unspecified parts of tongue': 11, 'Skin': 12,
			  'Corpus uteri': 13, 'Blood': 14, '0': 15, 'Cervix uteri': 16, 'Bladder': 17,
			  'Retroperitoneum and peritoneum': 18, 'Larynx': 19, 'Stomach': 20, 'Pancreas': 21,
			  'Bronchus and lung': 22, 'Thyroid gland': 23, 'Tonsil': 24, 'Rectum': 25,
			  'Hypopharynx': 26, 'Hematopoietic and reticuloendothelial systems': 27,
			  'Heart, mediastinum, and pleura': 28, 'Ovary': 29, 'Uterus, NOS': 30,
			  'Other and ill-defined sites in lip, oral cavity and pharynx': 31,
			  'Eye and adnexa': 32, 'Liver and intrahepatic bile ducts': 33,
			  'Connective, subcutaneous and other soft tissues': 34, 'Testis': 35,
			  'Rectosigmoid junction': 36, 'Other and unspecified parts of mouth': 37,
			  'Other and ill-defined sites': 38, 'Floor of mouth': 39,
			  'Other endocrine glands and related structures': 40, 'Lip': 41,
			  'Small intestine': 42, 'Gum': 43, 'Oropharynx': 44,
			  'Peripheral nerves and autonomic nervous system': 45,
			  'Bones, joints and articular cartilage of other and unspecified sites': 46,
			  'Unknown primary site': 47, 'Palate': 48,
			  'Bones, joints and articular cartilage of limbs': 49,
			  'Other and unspecified parts of biliary tract': 50, 'Meninges': 51,
			  'Gallbladder': 52,
			  'Spinal cord, cranial nerves, and other parts of central nervous system': 53,
			  'Other and unspecified male genital organs': 54}
	df = df[df['label'].map(df['label'].value_counts()) >= 600]

	df.label = [labels[item] for item in df.label]
	df['label'].replace('', np.nan, inplace=True)
	df.dropna()
	tumor_count = df.loc[df.label != 1].shape[0]
	normal_count = df.loc[df.label == 0].shape[0]
	logger.info("{} Normal samples, {} Tumor samples ".format(normal_count,tumor_count))
	columns = ['file_id','label']
	return df[columns]

if __name__ == '__main__':


	data_dir ="/Users/xin/Downloads/EE542/Labs/lab10/data/"
	# Input directory and label file. The directory that holds the data. Modify this when use.
	dirname = data_dir + "miRNA"
	label_file = data_dir + "files_meta.tsv"
	case_file = data_dir + "cases_meta.tsv"
	#output file
	outputfile = data_dir + "miRNA_matrix.csv"

	# extract data
	matrix_df = extractMatrix(dirname)
	label_df = extractLabel(label_file, case_file)

	#merge the two based on the file_id
	result = pd.merge(matrix_df, label_df, on='file_id', how="left")
	#print(result)

	#save data
	result.to_csv(outputfile, index=False)
	#print (labeldf)

 




