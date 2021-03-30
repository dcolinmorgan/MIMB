import pandas as pd
blood=pd.read_csv('~/analyses/LTRC/blood.txt', sep=" ")
lung=pd.read_csv('~/analyses/LTRC/lung.txt', sep=" ")

meta=pd.read_csv('~/analyses/LTRC/pheno_meta.txt', sep=" ")
copd=pd.read_csv('~/analyses/LTRC/modCopd_pathConservRna.pheno_022221.csv', sep=",")
ipf=pd.read_csv('~/analyses/LTRC/ipfRna.pheno_022221.csv', sep=",")
meta=meta.drop_duplicates()


jj=(pd.merge(pd.DataFrame(blood.columns),meta,left_on=0,right_on='new_id')).drop_duplicates()
hh=jj.drop_duplicates(subset=[0])[['ALIAS','new_id']]
hh.index=hh['ALIAS']
del hh['ALIAS']
hh=hh.T
hh.index.delete
dct={str(v[0]):k for k,v in hh.to_dict(orient='list').items()}
blood.rename(columns=dct, inplace=True)
# blood

jj=(pd.merge(pd.DataFrame(lung.columns),meta,left_on=0,right_on='new_id')).drop_duplicates()
hh=jj.drop_duplicates(subset=[0])[['ALIAS','new_id']]
hh.index=hh['ALIAS']
del hh['ALIAS']
hh=hh.T
hh.index.delete
dct={str(v[0]):k for k,v in hh.to_dict(orient='list').items()}
lung.rename(columns=dct, inplace=True)
# lung

cc=lung.corrwith(blood,drop=True,axis=0)
dd=lung.corrwith(blood,drop=True,axis=1)

cc.to_csv('~/analyses/LTRC/subject_corr.txt',sep='\t')
dd.to_csv('~/analyses/LTRC/cg_corr.txt',sep='\t')



blood_copd=blood[set(blood.columns) & set(copd[copd['modCopd_pathConserv']==1]['patid'])]
blood_ipf=blood[set(blood.columns) & set(ipf[ipf['ipf.clinpath']==1]['patid'])]
blood_control=blood[set(blood.columns) & set(ipf[ipf['ipf.clinpath']==0]['patid'])]

# aa=lung.corrwith(blood_copd,drop=True,axis=1)
bb=lung.corrwith(blood_ipf,drop=True,axis=1)
ccc=lung.corrwith(blood_control,drop=True,axis=1)
# aa.to_csv('~/analyses/LTRC/copd_cg_corr.txt',sep='\t')
bb.to_csv('~/analyses/LTRC/ipf_cg_corr.txt',sep='\t')
ccc.to_csv('~/analyses/LTRC/control_cg_corr.txt',sep='\t')



aa=lung.corrwith(blood_copd,drop=True,axis=0)
bb=lung.corrwith(blood_ipf,drop=True,axis=0)
ccc=lung.corrwith(blood_control,drop=True,axis=0)
aa.to_csv('~/analyses/LTRC/copd_subj_corr.txt',sep='\t')
bb.to_csv('~/analyses/LTRC/ipf_subj_corr.txt',sep='\t')
ccc.to_csv('~/analyses/LTRC/control_subj_corr.txt',sep='\t')
