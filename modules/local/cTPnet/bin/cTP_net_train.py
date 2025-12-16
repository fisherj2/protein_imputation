
print('loading packages')
import numpy as np
import torch.nn as nn
import pandas as pd
import torch.nn.functional as F
import torch
import torch.optim as optim
import sys
from scipy import stats
import argparse
import os
import gc
print('loading args')

#set seed
seed = 123
torch.manual_seed(seed)
torch.cuda.manual_seed(seed) 
torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU
np.random.seed(seed)


interactive=False

parser = argparse.ArgumentParser( prog = 'Script to train scipenn on scRNAseq data')
parser.add_argument('-d', '--basedir', required=True, help="pipeline base directory")
parser.add_argument('-b','--bench',  help='<Required> Set flag for benchmarking', required=True)
parser.add_argument('-f','--files', nargs='+', help='<Required> Set flag', required=True)

args = parser.parse_args()
dobenchmark = args.bench
input_files = args.files

#if length is one, probably need to split
if len(input_files)==1 or isinstance(input_files, (str)):
    print('adjusting input format')
    if isinstance(input_files, (list)):
        input_files=input_files[0]
    
    input_files=input_files.split('.csv')
    input_files=[s + '.csv' for s in input_files]

print(input_files)

#need to clean inputs of wrong character
characters_to_remove = ['[', ']', ',']
translation_table = str.maketrans('', '', ''.join(characters_to_remove))

rna_train_data_file = [file for file in input_files if 'training_data_rna_raw' in file][0]
rna_train_data_file = rna_train_data_file.translate(translation_table)

prot_train_data_file = [file for file in input_files if 'training_data_prot_raw' in file][0]
prot_train_data_file = prot_train_data_file.translate(translation_table)



rna_train_data_file=rna_train_data_file.strip()
prot_train_data_file=prot_train_data_file.strip()


print('files found')
#temp path for testing

#rna_train_data_file = '/fcrbiouatappn01/resbioinfo/data/MiroBio/shared/home/jfisher/analyses/protein_prediction/formal_analysis/output/training_files/training_data_rna_raw.csv'
#prot_train_data_file = '/fcrbiouatappn01/resbioinfo/data/MiroBio/shared/home/jfisher/analyses/protein_prediction/formal_analysis/output/training_files/training_data_prot_raw.csv'

X_list=[rna_train_data_file]
y_list=[prot_train_data_file]
header_list=['training_data']

#path to save model to
#loc=args.launchdir + '/output/cTPnet/'
# 
# if not os.path.exists(loc):
#     os.makedirs(loc)



#loading
X_final = pd.read_csv(X_list[0], index_col=0)
y_final = pd.read_csv(y_list[0], index_col=0)
protein_list=y_final.index.tolist()

del X_list
del y_list
gc.collect()

#need to clean inputs of wrong character
characters_to_remove = ['.']
translation_table = str.maketrans('', '', ''.join(characters_to_remove))


for i in range(len(protein_list)):
  val=protein_list[i]
  val=val.translate(translation_table)
  protein_list[i]=val


y_final.index=protein_list

repi=0

print('protein names formatted')

#
# gene_list=None
# X_final=None
# y_final=None
# for i,X_file in enumerate(X_list):
# 	print(X_file)
# 	X = pd.read_csv(X_file)
# 	if i==2:
# 		y = pd.read_csv(y_list[i],sep='\t')
# 	else:
# 		y = pd.read_csv(y_list[i])
# 	if interactive:
# 		X=X.transpose().sample(frac=0.1,random_state=4905).transpose()
# 	# Dealing with X's dimensionality
# 	if i==3 or i==4:
# 		gene=X.index.tolist()
# 		gene=[x[0:(len(x)-2)] for x in gene]
# 		X.index=gene
# 		gene=set(gene)
# 	else:
# 		gene = set(X.index.tolist())
# 	if gene_list is None:
# 		gene_list=gene
# 	else:
# 		gene_list=set(gene_list).intersection(gene)
# 	gene_list=list(gene_list)
# 	gene_list.sort()
# 	X=X.loc[gene_list,]
# 	if not X_final is None:
# 		X_final=X_final.loc[gene_list,]
# 	# Dealing with Y's dimensionality
# 	if i==2:
# 		protein=y.index.tolist()
# 		protein[protein.index('CD8a')]='CD8'
# 		protein[protein.index('CD127-IL7Ra')]='CD127'
# 		protein[protein.index('HLA.DR')]='HLA-DR'
# 		protein[protein.index('CD197-CCR7')]='CD197'
# 		protein[protein.index('CD278-ICOS')]='CD278'
# 	elif i==3 or i==4:
# 		protein=y.index.tolist()
# 		protein=[x.split("__")[0] for x in protein]
# 	else:
# 		protein = y.iloc[:,0].values.tolist()
# 		y=y.drop(columns=y.columns[0])
# 		if i==1:
# 			protein[protein.index('CCR5')]='CD195'
# 			protein[protein.index('CCR7')]='CD197'
# 	print(protein)
# 	y.index=protein
# 	y=y[X.columns]
# 	y=y.loc[protein_list,]
# 	# Add header to cell
# 	X.columns=list(map(lambda x: header_list[i]+'-'+x, X.columns.tolist()))
# 	y.columns=list(map(lambda x: header_list[i]+'-'+x, y.columns.tolist()))
# 	if i==0:
# 		X_final=X
# 		y_final=y
# 	else:
# 		X_final=pd.concat([X_final,X], axis=1)
# 		y_final=pd.concat([y_final,y], axis=1)


# Normalize y
shared_proteins=y_final.apply(lambda x: not x.isna().any(),axis=1)
y=y_final.apply(lambda x: np.log((x+1.0)/stats.gmean(x[shared_proteins]+1.0)), axis=0)
del(y_final)
# Use direct UMI count is hard to transfer across experiment
# Let's try normalize the data with seurat like method
X=X_final.apply(lambda x: np.log((x*10000.0/sum(x))+1.0), axis=0)
del(X_final)
# random cell order
X=X.T
X=X.sample(frac=1,random_state=4905)
# separate test data from train data
if interactive:
	X_test=X.sample(n=918,random_state=4905)# Need change after test
else:
	X_test=X.sample(n= round(X.shape[0] *0.1),random_state=4905)# Need change after test

y_test=y[X_test.index]

X=X.drop(X_test.index)
y=y.drop(columns=y_test.columns)
y=y[X.index]


# covert to tenor
X=torch.tensor(X.values)
X=X.type(torch.FloatTensor)
y=torch.tensor(y.values)
y=y.type(torch.FloatTensor)
y=torch.t(y)

X_test=torch.tensor(X_test.values)
X_test=X_test.type(torch.FloatTensor)
y_test=torch.tensor(y_test.values)
y_test=y_test.type(torch.FloatTensor)
y_test=torch.t(y_test)
n_batches=32


class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        #self.fc1 = nn.Linear(14505, 1000)
        self.fc1 = nn.Linear(X.shape[1], 1000)
        self.fc2 = nn.Linear(1000, 256)
        self.fc3 = nn.ModuleDict({})
        for p in protein_list:
        	self.fc3[p]=nn.Linear(256, 64)
        self.fc4 = nn.ModuleDict({})
        for p in protein_list:
        	self.fc4[p]=nn.Linear(64, 1)
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        outputs={}
        for p in protein_list:
        	outputs[p]=self.fc4[p](F.relu(self.fc3[p](x)))
        return outputs

net = Net()
# if repi==4:
# 	net.load_state_dict(torch.load('../model_optimal/training_04262020model_rep'+str(repi)+'_ep29'))
# else:
# 	net.load_state_dict(torch.load('../model_optimal/training_04262020model_rep'+str(repi)+'_ep54'))


criterion = nn.MSELoss()
optimizer = optim.Adam(net.parameters(), lr=0.001,amsgrad=True, weight_decay=0.001)
# optimizer = optim.Adagrad(net.parameters(), lr=lr_vec[repi], lr_decay=0.001)

max_epochs=200

train_loss=pd.DataFrame(np.zeros(shape=(len(protein_list),max_epochs)),index=protein_list)
test_loss=pd.DataFrame(np.zeros(shape=(len(protein_list),max_epochs)),index=protein_list)

# Init early stop
patience=30
best_score=None
Dy=len(protein_list)
estop_counter=pd.Series(np.zeros(Dy),index=protein_list)
early_stop=pd.Series([False]*Dy,index=protein_list)

for epoch in range(max_epochs):
	if all(early_stop):
		break
	running_loss=pd.Series(np.zeros(Dy),index=protein_list)
	for i in range(int(y.shape[0]/n_batches)):
		# Local batches and labels
		local_X, local_y = X[i*n_batches:min((i+1)*n_batches,X.shape[0]-1),], y[i*n_batches:min((i+1)*n_batches,y.shape[0]-1),]
		# zero the parameter gradients
		optimizer.zero_grad()
		# forward + backward + optimize
		outputs_dict = net(local_X)
		loss=None
		loss_count=0.0
		for p in protein_list:
			notNaN=(local_y[:,protein_list.index(p):(protein_list.index(p)+1)]==local_y[:,protein_list.index(p):(protein_list.index(p)+1)])
			loss_p=criterion(outputs_dict[p][notNaN],local_y[:,protein_list.index(p):(protein_list.index(p)+1)][notNaN])
			if not torch.isnan(loss_p):
				loss_count+=1.0
				running_loss[p]+=loss_p.item()
				if loss is None:
					loss=loss_p
				else:
					loss=loss+loss_p
		loss.backward()
		optimizer.step()
		if(i==(int(y.shape[0]/n_batches)-1)):
			train_loss.iloc[:,epoch]=(running_loss / 150)
		if i % 150 == 149:    # print every mini-batches
			print('[%d, %5d] loss: %.3f' % (epoch + 1, i + 1, sum(running_loss / 150)))
			running_loss=pd.Series(np.zeros(Dy),index=protein_list)
			sys.stdout.flush()
	test_outputs = net(X_test)
	test_outputs = [test_outputs[p] for p in protein_list]
	test_outputs=torch.transpose(torch.stack(test_outputs),0,1).view(X_test.shape[0],-1)
	test_loss_i=pd.Series([criterion(test_outputs[:,pi][y_test[:,pi]==y_test[:,pi]], y_test[:,pi][y_test[:,pi]==y_test[:,pi]]).item() for pi in range(Dy)],index=protein_list)
	test_loss.iloc[:,epoch]=test_loss_i
	#if epoch % 10 == 9:
		# f,ax=plt.subplots(figsize=(6,6))
		# ax.scatter(y_test.detach().numpy(),test_outputs.detach().numpy())
		# ax.plot([-2,5],[-2,5],ls='--',c='.3')
		# #ax.text(3,-2,'correlation: '+str(np.corrcoef(test_outputs.detach().numpy().flatten(),y_test.detach().numpy().flatten())[1,0]))
		# df = pd.DataFrame({"y_pred":test_outputs.detach().numpy().flatten(),'y_truth':y_test.detach().numpy().flatten()})
		# ax.text(3,-2,'correlation: '+str(round(df.corr().values[1,0],4)))
		# fig = ax.get_figure()
		# fig.savefig(loc+'figure_rep'+str(repi)+'_ep'+str(epoch)+'.pdf')
		# sys.stdout.flush()
		# plt.close(fig)
	if epoch % 5 == 4:
		torch.save(net.state_dict(), 'model_rep'+str(repi)+'_ep'+str(epoch))
	# Implement early stopping
	if best_score is None:
		best_score=test_loss_i
	else:
		for p in protein_list:
			if test_loss_i[p]>(best_score[p]-0.001) and (not early_stop[p]):
				estop_counter[p]+=1
				if estop_counter[p]>=patience:
					early_stop[p]=True
			else:
				best_score[p]=test_loss_i[p]
				estop_counter[p]=0
		print(estop_counter)

print('Finished Training')

#save to directory
#torch.save(net.state_dict(), loc+'final_model_rep'+str(repi)+'_ep'+str(epoch))
train_loss.index=['train_'+p for p in protein_list]
test_loss.index=['test_'+p for p in protein_list]
log=pd.concat([train_loss,test_loss])



#save to working directory, pass forward along pipe



if dobenchmark=='true':
    #get benchmarking prefix

	basename=os.path.basename(rna_train_data_file)
	prefix= basename.replace("_training_data_rna_norm.csv", "") 
	#torch.save(net.state_dict(), loc+prefix +'final_model_rep'+str(repi)+'_ep'+str(epoch))
	torch.save(net.state_dict(), prefix+ 'final_model_rep'+str(repi)+'_ep'+str(epoch))
	#log.to_csv(loc+prefix+'log_rep'+str(repi)+'.csv')
else:
    #save for later
    #torch.save(net.state_dict(), loc+'final_model_rep'+str(repi)+'_ep'+str(epoch))
    torch.save(net.state_dict(), 'final_model_rep'+str(repi)+'_ep'+str(epoch))
    #log.to_csv(loc+'log_rep'+str(repi)+'.csv')
