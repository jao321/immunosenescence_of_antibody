clonotyped = glob.glob("path/to/folder/with/tsv_YClon_clonotyped_files/**/*YClon_clonotyped.tsv", recursive=True)

ID_lst =[] 
matu =[]
for x in clonotyped:
    print(os.path.basename(x))
    df = pd.read_csv(x, sep="\t")
    clone_id_group = df.groupby(['clone_id'])
    result = pd.DataFrame()
    result["number_of_v_gene"] = clone_id_group.v_call.nunique()
    result["number_of_j_gene"] = clone_id_group.j_call.nunique()
    result["ID"] = os.path.basename(x)
    tmp = pd.DataFrame()
    tmp["v_j_product"] = clone_id_group.j_identity.apply(lambda x: x*0.2) + clone_id_group.v_identity.apply(lambda x: x*0.8)
    tmp = tmp.reset_index()
    result["v_j_product"] = tmp.groupby("clone_id").v_j_product.max()
    result = result.reset_index()
    ID_lst.append(os.path.basename(x).split("_")[0])
    matu.append(len(result[result["v_j_product"]>=0.9])/len(df["clone_id"].unique()))
maturity = pd.DataFrame({'ID':ID_lst,'matu':matu})