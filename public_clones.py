import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import upsetplot

public = pd.DataFrame()
with pd.read_csv("path/to/file/tsv_YPub_input_YClon_clonotyped.tsv",
                 sep="\t",
                 chunksize=1000000) as public:
    for chunk in public:
        public = pd.concat([public,chunk[['sequence_id','productive', 'v_call',
                                           'origin_repertoire', 'clone_id','clone_seq_count']]])


meta = pd.read_csv("path/to/file/with/metadata.csv")

public = pd.merge(public_covidao,meta,on='sample_ID')

Control_set=set(public[((public.CITY=="BH") & (public['CLINICAL CLASSIFICATION  ']=="Control"))].clone_id.unique())
Mild_set=set(public[((public.CITY=="BH") & (public['CLINICAL CLASSIFICATION  ']=="Mild"))].clone_id.unique())
Hospitalized_set=set(public[((public.CITY=="BH") & (public['CLINICAL CLASSIFICATION  ']=="Hospitalized"))].clone_id.unique())
GV_set=set(public[((public.CITY=="GV") & (public['CLINICAL CLASSIFICATION  ']=="Mild"))].clone_id.unique())

from matplotlib.ticker import FuncFormatter
fig = plt.figure(figsize=(20, 10))
set_names = ['Control','Mild_EA','Mild-NEA','Hospitalized_NEA']
all_clones = Control_set.union(GV_set).union(Mild_set).union(Hospitalized_set)
df = pd.DataFrame([[e in Control_set, e in GV_set, e in Mild_set, e in Hospitalized_set] for e in all_clones], columns = set_names)
df_up = df.groupby(set_names).size()
upsetplot.plot(df_up,orientation='horizontal', fig=fig, element_size=None,show_counts='{:,d}')
ax = plt.gca()
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: '{:,}'.format(int(x))))

plt.rcParams.update({'font.size':12})

plt.savefig("path/to/samples_analysis/public_clones.png")
