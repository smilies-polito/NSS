#
# Copyright © 2022 Politecnico di Torino, Control and Computer Engineering Department, SMILIES group
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
# associated documentation files (the "Software"), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial 
# portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
import os.path

import sys
import os
import time

import numpy as np
import pandas as pd
import seaborn as sns

import plotly.express as px
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import gridspec
from sklearn.model_selection import train_test_split
from sklearn.inspection import permutation_importance
from scipy.cluster.hierarchy import dendrogram
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering

pd.options.mode.chained_assignment = None  # default='warn'


#features from literature which are deemed important for neuron AP analysis
features_liter = ['tau', 'vrest', 'threshold_v_ramp', 'upstroke_downstroke_ratio_ramp', 'peak_v_ramp', 'fast_trough_v_ramp']


def main(dataset_type='PatchSeqDataset', cell_type='aspiny',  stimulation='ramp', n_clust_to_search=2, vis_label='CLUSTER'):
    """
    Main Function of the script. It performs:
    - Data loading (either PatchClampDataset or PatchSeqDataset)
    - Data preprocessing
        --Feature selection and extraction
        --Feature inspection (histograms, boxplots)
    - Reduction and visualization through PCA according to different criteria (neuron family type, clusterization)
    - Clusterization for the recognition of neurons subgroups using Agglomerative clustering
    - Inspection of features in clusters and their relevance to clusterization
    - Save line of neurons in clusters to correlate neuron type with ef features of the spike

    @param dataset_type: dataset to load --> 'PatchClampDataset' or 'PatchSeqDataset'
    @param cell_type: cells to select from the dataset --> 'full' (whole dataset),  'aspiny' (inhibitors), 'spiny' (excitatory), 'sparsely spiny' (mixed behaviour) (string)
    @param stimulation: type of electro-physiological stimulation to explore --> 'ramp', 'short_square' or 'long_square' (string)
    @param n_clust_to_search: user defined clusters to search for the data (default to 2, fix according to automatic number of cluster selection methods)-->  2 <= 'n_clust_to_search' < +inf (integer)
    @param vis_label: correlate analysed features according to cell_type (inhibitor/excitatory) or identified clusters -->'NEURON_TYPE' or 'CLUSTER' (string)
    """

    global path_for_images, path_for_csv, path_for_misc, clusters
    

    #paths to save images and files generated
    path_for_images=os.path.join('..', 'output', dataset_type, 'IMAGES')
    path_for_csv=os.path.join('..', 'output', dataset_type)
    path_for_misc=os.path.join('..', 'output', dataset_type, 'MISC')
    
    
    if not os.path.isdir(path_for_images):
        os.makedirs(path_for_images)
        
    if not os.path.isdir(path_for_csv):
        os.makedirs(path_for_csv)
        
    if not os.path.isdir(path_for_misc):
        os.makedirs(path_for_misc)

    
    
    clusters=n_clust_to_search

    stimulation="_"+stimulation

    # select features from PatchClampDataset dataset according to required 'stimulation'
    
    for i in range(0, len(features_liter)):
        if '_ramp' in features_liter[i]:
            features_liter[i]=features_liter[i].replace('_ramp', stimulation)

    # Load the data
    if dataset_type=='PatchClampDataset':
       dataset_neurons=load_dataPatchClampDataset(cell_type)
    else:
       dataset_neurons = load_dataPatchSeqDataset(cell_type)

    # Preprocessing
    dataset_neurons_clean, neuron_type, donor, region, specimen, transcript_ids, cell_type_ids = preprocessData(dataset_neurons, dataset_type=dataset_type,
                                                                                type_of_input=stimulation)

    # Map categorical properties to numbers for grouping/visualization purposes
    neuron_type = neuron_type.map({"aspiny": 0, "spiny": 1, "sparsely spiny": 2})

    # Compute PCA components of the preprocessed Dataset for visualization
    PCs = computePCA(dataset_neurons_clean)


    #Rename
    PCs = PCs[['pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6', 'pc7']]

    # OPT1: silhoutte score to identify optimal number of clusters k
    CLUSTERING_METHOD='Agglomerative'

    optimal_k=silhoutte_plot(dataset_neurons_clean, CLUSTERING_METHOD)

    if n_clust_to_search==-1:
        n_clust_to_search=optimal_k


    if CLUSTERING_METHOD=='Agglomerative':
        #dendogram_construction(dataset_neurons_clean)
        model=AgglomerativeClustering(n_clusters=n_clust_to_search).fit(dataset_neurons_clean)


    # Visualization with respect to agglomerative clusterization (user-defined k)
    clusterization = model
    if(n_clust_to_search==3):
        if dataset_type=='PatchClampDataset':
            cmap={'EC0':"#440154",'EC1':"#21908c", 'EC2':"#fde725"}   #1,2,0
            visualizationPCA(PCs, list(map({1:'EC0', 2:'EC1', 0:'EC2'}.get, clusterization.labels_.tolist(),)), 3, cmap, cell_labels='NSS_labels')
        else:
            cmap = {'EC0': "#440154", 'EC2': "#fde725", 'EC1': "#21908c"}  # 1,2,0
            visualizationPCA(PCs, list(map({1: 'EC0', 2: 'EC2', 0: 'EC1'}.get, clusterization.labels_.tolist(), )), 3, cmap, cell_labels='NSS_labels')
    else:
        if dataset_type=='PatchClampDataset':
            cmap={'EC0':"#440154",'EC1':"#fde725"}
            visualizationPCA(PCs, list(map({1:'EC0', 0:'EC1'}.get, clusterization.labels_.tolist())), 3, cmap, cell_labels='NSS_labels')
        else:
            cmap = {'EC0': "#440154", 'EC1': "#fde725"}  # 1,2,0
            visualizationPCA(PCs, list(map({1: 'EC0', 0: 'EC1'}.get, clusterization.labels_.tolist(), )), 3,
                             cmap, cell_labels='NSS_labels')

    selected_feat = list(dataset_neurons_clean.columns)

    # According to 'vis_label', correlate the obtained results to
    # clusterization or 'neuron_type'  (NB: neuron type only if 'cell_type' is set to 'full')

    if (vis_label == 'CLUSTER'):
       dataset_neurons_clean['LABEL'] = clusterization.labels_
    else:
       dataset_neurons_clean['LABEL'] = neuron_type


    if dataset_type=='PatchSeqDataset':

        #Save transcriptomics ids for members of clusters
        save_clusters_transcr(dataset_neurons_clean, transcript_ids, cell_type_ids)
    else:
        
        # Save specimen line to identify members of clusters
        save_clusters_line(dataset_neurons_clean, specimen, cell_type_ids)
        
        # saving NSS clusters 
        save_clusters_transcr(dataset_neurons_clean, specimen, cell_type_ids)

        
        directory = os.path.join(path_for_csv,"cre_lines_all")
        
        if not os.path.isdir(directory):
            os.makedirs(directory)
        
        with open(os.path.join(path_for_csv,"cre_lines_all", "cre_lines.csv"), 'w') as file:
            for s, id in zip(specimen, cell_type_ids):
                file.write(str(s) + "," + str(id) + "\n")

        # create the cre lines - cell types association dictionary
        tot = []
        
        directory = os.path.join(path_for_csv, "cre_lines_clusters")
        
        if not os.path.isdir(directory):
            os.makedirs(directory)
        
        directory = os.path.join(path_for_csv, "cre_lines_clusters")
        for filename in os.listdir(directory):
            f = os.path.join(directory, filename)
            f_out = os.path.join(path_for_misc, filename.rstrip(".csv") + "_type.csv")
            list_pop = cluster_population_lines(f)
            tot.append(list_pop)
            # translate the cre lines into cell types
            population_change_name(f, f_out)

        population_change_name(os.path.join(path_for_csv,"cre_lines_all", "cre_lines.csv"), os.path.join(path_for_misc, "cre_lines_type.csv"))
        # percentage population for clusters
        for i in range(len(tot[0]) - 1):
            pop = 0
            print("Population", i, " is divided in clusters as:")
            for j in range(len(tot)):
                pop = pop + tot[j][i]
            for j in range(len(tot)):
                print(tot[j][i] / pop * 100, "% in cluster", j)


    # Visualization with respect to neuron type
    if dataset_type=="PatchClampDataset":

        with open(os.path.join(path_for_misc, 'cre_lines_type.csv'), 'r') as file:
            cre_lines=[l.strip('\n').split(",")[0] for l in file.readlines()]

        cre_lines=pd.Series(cre_lines)

        #Visualization with respect to specimen only for PatchClampDataset
        visualizationPCA(PCs, cre_lines, 3, {'Pvalb':'#53B400', 'Sst':'#A58AFF', 'Vip/Lamp5':'#FB61D7', 'Unsure':'#D3D3D3'}, cell_labels='Cre_lines_labels')

    # Violin plots of features according to either clusterization or 'neuron type'  (NB: neuron type only if 'cell_type' is set to 'full')
    #violin_plots(dataset_neurons_clean, selected_feat, vis_label)

    # Correlation to clusterization or neuron type
    # (we use Spearman since most of features appear non-normally distributed from violin and histograms)
    p_correlation = dataset_neurons_clean.corr(method="spearman")

    sns.set(font_scale=1.5)
    upper_tri = (p_correlation).where(np.triu(np.ones(p_correlation.shape), k=1).astype(bool))
    upper_tri = upper_tri.iloc[:-1, 1:]

    # Create custom layout using GridSpec
    fig = plt.figure(figsize=(10, 8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[15, 1], height_ratios=[1, 15])
    ax_heatmap = plt.subplot(gs[1, 0])
    ax_cbar = plt.subplot(gs[0, 0])

    # Plot correlation matrix
    heatmap = sns.heatmap(upper_tri, ax=ax_heatmap, xticklabels=True, yticklabels=True, vmin=-1, vmax=1, annot=True,
                          cmap=sns.color_palette("icefire", as_cmap=True), cbar=False)


    # Add colorbar
    cbar = fig.colorbar(heatmap.get_children()[0], cax=ax_cbar, orientation='horizontal', pad=0.02)
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    cbar.set_label('Correlation')

    # Rotate x-axis tick labels for better readability
    ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels(), rotation=70)
    
        
    plt.tight_layout()
    plt.savefig(os.path.join(path_for_images, 'correlation_matrix_k'+str(clusters)+'.png'), dpi=300, bbox_inches='tight')
    #plt.show()

    # Features importance using random forest
    y = dataset_neurons_clean['LABEL']
    X = pd.DataFrame.to_numpy(dataset_neurons_clean.drop(['LABEL'], axis='columns', inplace=False))
    RF_features_importance(X, y, selected_feat)

    # show distribution of features in the identified clusters
    if n_clust_to_search==3:
        cmap=["#440154", "#21908c", "#f2c705"]
        legend=["EC2", "EC1", "EC0"]
    elif n_clust_to_search==2:
        cmap=["#440154", "#f2c705"]
        legend=['EC1', 'EC0']

    clusters_inspection(dataset_neurons_clean, selected_feat, cmap, legend)


def save_clusters_transcr(dataset, transcr_ids, celltype_ids):
    """
        Function to output the transcriptomics ids of the neurons in the dataset grouped according to LABEL feature
        NB: LABEL can be either cell_type (aspiny, spiny, sparsely spiny) or label of automatic clusterization (0,1,2..).
        The output are stored in file "transcriptomics_clusters.csv" in "output" folder

        @param dataset: dataset containing the neurons (Pandas Dataframe)
        @param transcr_ids: transcriptomic ids of the neurons in the dataset (Pandas Series)
        """

    if not os.path.isdir(os.path.join(path_for_csv, "NSS_clusters")):
            os.makedirs(os.path.join(path_for_csv, "NSS_clusters"))
            
    with open(os.path.join(path_for_csv, "NSS_clusters", "NSS_clusters_k"+str(clusters)+".csv"), 'w') as file:
        
        
        for l in dataset['LABEL'].unique():

            # save specimen line to identify members of clusters
            transcr_cl = transcr_ids[dataset['LABEL'] == l]

            dataset_ridotto=dataset[dataset['LABEL']==l]
            dataset_ridotto.drop(['LABEL'], axis='columns', inplace=True)

            #I look for cluster median cells
            median_cell=dataset_ridotto.median(axis=0, numeric_only=True)
            min_dist = (-1, 0)
            for (index,row) in dataset_ridotto.iterrows():
                dist = (np.linalg.norm(row - median_cell))
                if(min_dist[0]==-1 or min_dist[1]>dist):
                    min_dist=(index, dist)

        

            median_cell_id=celltype_ids[int(min_dist[0])]
            print("Median cell of NSS cluster "+str(l)+": "+str(median_cell_id))

            for s in transcr_cl:
                file.write(s + ","+str(l)+'\n')


def save_clusters_line(dataset, specimen, celltype_ids):
    """
    Function to output the family line of the neurons in the dataset grouped according to LABEL feature
    NB: LABEL can be either cell_type (aspiny, spiny, sparsely spiny) or label of automatic clusterization (0,1,2..).
    The output are stored in .txt files named using the notation "line_cl"+LABEL_VALUE+".csv"

    @param dataset: dataset containing the neurons (Pandas Dataframe)
    @param specimen: family line of the neurons in the dataset (Pandas Series)
    """
    for l in dataset['LABEL'].unique():
        # save specimen line to identify members of clusters
        specimen_cl = specimen[dataset['LABEL'] == l]
        celltype_ids_cl=celltype_ids[dataset['LABEL']==l]
        dataset_ridotto = dataset[dataset['LABEL'] == l]
        dataset_ridotto.drop(['LABEL'], axis='columns', inplace=True)

        # I look for cluster median cells
        median_cell = dataset_ridotto.median(axis=0, numeric_only=True)
        min_dist = (-1, 0)
        for (index, row) in dataset_ridotto.iterrows():
            dist = (np.linalg.norm(row - median_cell))
            if (min_dist[0] == -1 or min_dist[1] > dist):
                min_dist = (index, dist)
        
        median_cell_id = celltype_ids_cl[int(min_dist[0])]
        print("Median cell of cre lines cluster " + str(l) + ": " + str(median_cell_id))

        directory = os.path.join(path_for_csv, "cre_lines_clusters")
        
        if not os.path.isdir(directory):
            os.makedirs(directory)
            
        with open(os.path.join(directory, "lines_k"+str(clusters)+"_cl"+str(l)+".csv"), 'w') as file:
            for s, id  in zip(specimen_cl, celltype_ids_cl):
                file.write(str(s)+","+str(id)+ "\n")
                
        directory = os.path.join(path_for_csv, "NSS_clusters")
        if not os.path.isdir(directory):
            os.makedirs(directory)

def violin_plots(dataset, selected_feat, vis_label):
    """
    Function to visualize distribution of selected features of dataset using violin plots. Current version assumes
    the existence of two groups to analyse separately according to feature 'LABEL' (either 0 or 1).

    @param dataset: dataset of neurons to analyse (Pandas Dataframe)
    @param selected_feat: list of features in dataset to explore through violin plots (NB: pass a list also for single feature - e.g., ['deep_potential']_) (list of string)
    @param vis_label: string containing meaning of 'LABEL' property in dataset (either 'cell_type' or label of automatic clusterization) (string)
    """
    df_melt = pd.melt(dataset, id_vars=['LABEL'], value_vars=selected_feat, var_name='Features',
                      value_name='Normalized value')

    df_melt['LABEL'] = df_melt['LABEL'].apply(lambda x: 0 if x < 1 else 1)
    df_melt['LABEL'] = df_melt['LABEL'].map({0: "0", 1: "1"})
    plt.figure()
    plt.title("Violin plots discriminating by " + vis_label)
    my_pal = {"0": '#FFE135', "1": '#863B8D'}
    sns.set(rc={'font.weight': 'bold'})
    sns.set_style("whitegrid")
    ax = sns.violinplot(x='Features', y='Normalized value', hue="LABEL", data=df_melt, palette=my_pal, split=True,
                        kind='kde')
    ax.set_xticklabels(selected_feat, fontsize=13, rotation=45)
    plt.legend(fontsize='x-large', title_fontsize='40')
    #plt.savefig(os.path.join(path_for_images, 'violin_plots_', selected_feat, '.png'), dpi=300.0, bbox_inches='tight')
    #plt.show()

def clusters_inspection(dataset, feat_to_inspect=[], cmap=["#440154", "#21908c", "#f2c705"], labels=["EC2", "EC1", "EC0"]):
    """
    Function to visualize distribution of the features in the groups defined by property 'LABEL' through KDE plots
    @param dataset: dataset of neurons to inspect (Pandas Dataframe)
    @param feat_to_inspect: features to inspect (list of strings, also for single feature)
    """
    dataset_copy=dataset.copy(deep=True)
    
    #to have same colors as in scatter plot
    dataset_copy['LABEL']=dataset_copy['LABEL'].replace({0:1, 1:0})
    
    print("KDE plot generation...")
    for f in feat_to_inspect:
        print("Inspecting feature "+f)
        plt.figure()
        sns.kdeplot(data=dataset_copy, x=f, hue='LABEL', palette=cmap, linewidth=2.5)
        plt.legend(labels=labels)
        plt.savefig(os.path.join(path_for_images, 'KDE_plot_k'+str(clusters)+'_'+f), dpi=300.0, bbox_inches='tight')
        #plt.show()

def RF_features_importance(X, y, selected_feat):
    """
    Function to compute features importance for clusters using Random Forest classifier from
    Sklearn <https://scikit-learn.org/stable/auto_examples/ensemble/plot_forest_importances.html>

    @param X: matrix containing data to analyse, rows are datapoints, columns are features (Numpy matrix)
    @param y: label to infer for data in x (Numpy array)
    @param selected_feat: names of features to evaluate (list of strings)
    """
    
    print("Computing feature importance...")
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=42)

    # A random forest classifier will be fitted to compute the feature importances.
    feature_names = selected_feat
    forest = RandomForestClassifier(random_state=0)
    forest.fit(X_train, y_train)

    # %%
    # Feature importance based on mean decrease in impurity
    # -----------------------------------------------------
    # Feature importances are provided by the fitted attribute
    # `feature_importances_` and they are computed as the mean and standard
    # deviation of accumulation of the impurity decrease within each tree.
    #
    # .. warning::
    #     Impurity-based feature importances can be misleading for **high
    #     cardinality** features (many unique values). See
    #     :ref:`permutation_importance` as an alternative below.
    #start_time = time.time()
    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_], axis=0)
    #elapsed_time = time.time() - start_time

    #print(f"Elapsed time to compute the importances: {elapsed_time:.3f} seconds")


    # Let's plot the impurity-based importance.
    forest_importances = pd.Series(importances, index=feature_names)

    fig, ax = plt.subplots()
    forest_importances.plot.bar(yerr=std, ax=ax)
    ax.set_title("Feature importances using MDI")
    ax.set_ylabel("Mean decrease in impurity")
    fig.tight_layout()


    # Feature importance based on feature permutation
    # -----------------------------------------------
    # Permutation feature importance overcomes limitations of the impurity-based
    # feature importance: they do not have a bias toward high-cardinality features
    # and can be computed on a left-out test set.

    #start_time = time.time()
    result = permutation_importance(
        forest, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2
    )
    #elapsed_time = time.time() - start_time
    #print(f"Elapsed time to compute the importances: {elapsed_time:.3f} seconds")

    forest_importances = pd.Series(result.importances_mean, index=feature_names)


    # The computation for full permutation importance is more costly. Features are
    # shuffled n times and the model refitted to estimate the importance of it.
    # Please see :ref:`permutation_importance` for more details. We can now plot
    # the importance ranking.

    plt.rc('axes', titlesize=16)  # fontsize of the axes title
    plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=18)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)  # fontsize of the tick labels
    plt.rc('legend', fontsize=18)  # legend fontsize

    #save importance for plotting k=2 and k=3 cases together
    forest_importances.to_csv(os.path.join(path_for_misc, "RF_features_importance_k"+str(clusters)+".csv"))
    np.savetxt(os.path.join(path_for_misc, "RF_features_importance_std_k"+str(clusters)+".csv"), result.importances_std, delimiter=",");


    try:
    
        #it takes both analyses to be performed at least once to plot k=2 and k=3 together
        # so an exception is thrown if one of the files is not present
    
        if clusters==3:
            fi_to_load="RF_features_importance_k2.csv"
            fi_std_to_load="RF_features_importance_std_k2.csv"
        elif clusters==2 :
            fi_to_load = "RF_features_importance_k3.csv"
            fi_std_to_load = "RF_features_importance_std_k3.csv"
        fi_loaded=pd.read_csv(os.path.join(path_for_misc, fi_to_load));
        fi_loaded=fi_loaded.set_index(forest_importances.index)
        forest_importances_df=pd.concat([forest_importances*100, pd.Series(fi_loaded['0']*100, dtype='float32')], axis=1);
        forest_importances_df.columns=['K=2', 'K=3']

        fi_std_loaded=np.loadtxt(os.path.join(path_for_misc, fi_std_to_load))
        importances_std= np.append(np.reshape(result.importances_std*100, (10,1)), np.reshape(fi_std_loaded*100, (10,1)), axis=1)


        fig, ax = plt.subplots()
        forest_importances_df.plot.bar(yerr=np.transpose(importances_std), ax=ax)
        ax.set_title("Feature importances using permutation on full model", fontsize=15)
        ax.set_ylabel("Mean accuracy (%) decrease", fontsize=15)
        ax.set_yticks([0, 2, 4, 6, 8, 10, 12])
        fig.tight_layout()
    
    
        plt.savefig(os.path.join(path_for_images, 'bar_plots_mean_accuracy_decrease.png'), dpi=300.0, bbox_inches='tight')
        #plt.show()
    
    except FileNotFoundError:
        print("WARNING: Perform the chosen analysis with both k=2 and k=3 to generate the combined feature importance plot.")
    
    except:
        print("Something went wrong with the combined feature importance plot.")

    # %%
    # The same features are detected as most important using both methods. Although
    # the relative importances vary. As seen on the plots, MDI is less likely than
    # permutation importance to fully omit a feature.

def silhoutte_plot(dataset, clust_type='Agglomerative'):
    """
    Function to tune numbers_of_clusters for Agglomerative clustering using silhoutte score.
    @param dataset: dataset on which to compute the silhoutte score
    """
    
    print("Silhouette plot generation...")
    sil_scores = {}
    if clust_type=='Agglomerative':
        for k in range(2, 10):
            clusterization = AgglomerativeClustering(n_clusters=k).fit(dataset)
            sil_scores[k] = metrics.silhouette_score(dataset, clusterization.labels_)

    #print(sil_scores.keys())
    #print(sil_scores.values())

    plt.figure()
    plt.plot(list(sil_scores.keys()), list(sil_scores.values()))
    plt.xlabel("Number of clusters K")
    plt.ylabel("Silhouette score")
    
    
    optimal_k = max(sil_scores, key=lambda k: sil_scores[k])
    
    plt.savefig(os.path.join(path_for_images, 'silhouette.png'),dpi=300.0, bbox_inches='tight')
    #plt.show()
    return optimal_k


def plot_dendrogram(model, **kwargs):
    """
    Function to plot the dendogram of an Agglomerative Clustering

    @param model: agglomerative clustering model (sklearn.cluster.AgglomerativeClustering object)
    @param kwargs: parameters for dendogram renderings (list of objects)
    """
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)


def dendogram_construction(dataset):
    """
    Function to tune numbers_of_clusters for Agglomerative clustering using dendogram construction

    @param dataset: dataset on which to compute the dendogram (Pandas dataframe)
    """
    clusterization = AgglomerativeClustering(n_clusters=None, distance_threshold=0).fit(dataset)
    plt.title("Hierarchical Clustering Dendrogram")
    plot_dendrogram(clusterization, truncate_mode="level", p=3)
    plt.xlabel("Number of points in node (or index of point if no parenthesis).", fontsize=14)
    #plt.show()



def computePCA(dataset_n):
    """
    Function to compute Principal Components Analysis on a dataset. Variance explained bu each principal component is
    shown in figure

    @param dataset_n: dataset over which to compute PCA (already normalized using z-score) (Pandas dataframe)

    @return: PCs: original dataset transformed according to computed PCA (Pandas dataframe)
    """
    # Compute pca on the normalized data
    pca_cmp = PCA(len(dataset_n.columns))
    pca_cmp.fit(dataset_n)

    # Create the visualization plot of explained variance by each PC
    exp_var_pca = pca_cmp.explained_variance_ratio_
    cum_sum_eigenvalues = np.cumsum(exp_var_pca)

    #see explained variance by each PC
    # plt.bar(range(0, len(exp_var_pca)), exp_var_pca, alpha=0.5, align='center', label='Individual explained variance')
    # plt.step(range(0, len(cum_sum_eigenvalues)), cum_sum_eigenvalues, where='mid',
    #          label='Cumulative explained variance')
    # plt.ylabel('Explained variance ratio')
    # plt.xlabel('Principal component index')
    # plt.legend(loc='best')
    # plt.tight_layout()
    # plt.show()

    # Apply the transformation to original dataset
    PCs = pd.DataFrame(pca_cmp.transform(dataset_n))
    PCs.columns = ["pc" + str(i) for i in range(1, len(dataset_n.columns) + 1)]

    #visualize contributions of original features to PCs, normalized by PC importance (method from https://medium.com/@gaurav_bio/creating-visualizations-to-better-understand-your-data-and-models-part-1-a51e7e5af9c0)
    # feat_contr_to_pc=dataset_n.reset_index(drop=True).T.dot(PCs)
    # feat_contr_to_pc=(feat_contr_to_pc-feat_contr_to_pc.mean())/feat_contr_to_pc.std()
    # prova_vis=feat_contr_to_pc.multiply(exp_var_pca, axis=1)
    # plt.figure()
    # sns.heatmap(abs(prova_vis), xticklabels=True, yticklabels=True, annot=False, cmap=sns.cubehelix_palette(as_cmap=True))
    # plt.show()
    return PCs


def visualizationPCA(PCs, y, n_pc=2, cmap=None, cell_labels=''):
    """
    Function to visualize scatter of the data according to first two (or three) Principal Components

    @param PCs: dataframe containing Principal Components of the data (at least 2) (Pandas dataframe)
    @param y: label to apply to the data for coloring points (Numpy array or Pandas Series)
    @param n_pc: number of principal components to visualize (either 2 or 3) (integer)
    """
    
    label=[str(l) for l in y]

    if cmap!=None:
        if n_pc == 2:
            fig = px.scatter(PCs, x='pc1', y='pc2', color=label, color_discrete_map=cmap)
            fig.update_layout(font=dict(size=18))
            fig.write_html(os.path.join(path_for_images, str('NSS_embedding_' + cell_labels + '_k' + str(clusters) + '.html')))
            #fig.show()
        else:
            fig = px.scatter_3d(PCs, x='pc1', y='pc2', z='pc3', color=label, color_discrete_map=cmap)
            fig.update_layout(font=dict(size=18))
            fig.write_html(os.path.join(path_for_images, str('NSS_embedding_' + cell_labels + '_k' + str(clusters) + '.html')))
            #fig.show()
    else:

        if n_pc == 2:
            fig = px.scatter(PCs, x='pc1', y='pc2', color=label)
            fig.update_layout(font=dict(size=18))
            fig.write_html(os.path.join(path_for_images, str('NSS_embedding_' + cell_labels + '_k' + str(clusters) + '.html')))
            #fig.show()
        else:
            fig = px.scatter_3d(PCs, x='pc1', y='pc2', z='pc3', color=label)
            fig.update_layout(font=dict(size=18))
            fig.write_html(os.path.join(path_for_images, str('NSS_embedding_' + cell_labels + '_k' + str(clusters) + '.html')))
            #fig.show()


def preprocessData(dataset_neurons, dataset_type='PatchClampDataset', type_of_input='ramp'):
    """
    Function to preprocess data. It performs:
        -Features selection and extraction
        -Data cleaning (Nan and outlier removal)
        -preliminary data exploration (histograms and boxplots of features, correlation matrix)
        -data normalization through z-score

    @param dataset_neurons: input dataset of neurons to preprocess (Pandas dataframe)
    @param dataset_type: type of the dataset, either "PatchClampDataset" or "PatchSeqDataset" (string)
    @param type_of_input: specify shape of input stimulation to analyze (NB: for dataset_type=='PatchSeqDataset' only 'ramp' is valid) (string)

    @return dataset_neurons_znorm[features]: input dataset normalized, considering only selected features (Pandas Dataframe)
    @return type: type of the neurons included in preprocessed dataset (aspiny, spiny, sparsely spiny) (Series of strings)
    @return donor: donor of the neurons included in preprocessed dataset ('Human' or 'Mus musculus') (Pandas Series of strings)
    @return region: region of the brain of the neurons included in preprocessed dataset (Pandas Series of integer IDs)
    @return specimen: identifier of family line of the neurons included in preprocessed dataset (Pandas Series of integers)
    """
    
    if dataset_type == 'PatchClampDataset':
        print("PatchClampDataset pre-processing... ")

    if dataset_type == "PatchSeqDataset":
        if type_of_input!='_ramp':
            print("ERROR! For PatchSeqDataset only RAMP stimulation is available. Try again")
            exit(-1)
        print("PatchSeqDataset pre-processing... ")

    new_feat=10


    # add the derived features (indipendent from chosen dataset)
    dataset_neurons['AP_halfwidth'] = (dataset_neurons['peak_t'+type_of_input] - dataset_neurons['threshold_t'+type_of_input]) / 2
    dataset_neurons['Down_width'] = (dataset_neurons['fast_trough_t'+type_of_input] - dataset_neurons['peak_t'+type_of_input])
    dataset_neurons['UpDown_width'] = (dataset_neurons['fast_trough_t'+type_of_input] - dataset_neurons['threshold_t'+type_of_input])
    dataset_neurons['Width'] = abs((dataset_neurons['threshold_t'+type_of_input] + dataset_neurons['AP_halfwidth']) - (dataset_neurons['peak_t'+type_of_input] + ((dataset_neurons['peak_t'+type_of_input] - dataset_neurons['fast_trough_t'+type_of_input]) / 2)))
    dataset_neurons['Heigth'] = abs((dataset_neurons['peak_v'+type_of_input] - dataset_neurons['fast_trough_v'+type_of_input]))
    dataset_neurons['ΔV_deep'] = abs(dataset_neurons['threshold_v'+type_of_input] - dataset_neurons['fast_trough_v'+type_of_input])
    dataset_neurons['ΔV_ratio'] = (abs(dataset_neurons['peak_v'+type_of_input] - dataset_neurons['threshold_v'+type_of_input])) / abs(dataset_neurons['peak_v'+type_of_input] - dataset_neurons['fast_trough_v'+type_of_input])
    dataset_neurons['Slope_deep'] = (abs(dataset_neurons['fast_trough_v'+type_of_input] - dataset_neurons['threshold_v'+type_of_input])) / abs(dataset_neurons['fast_trough_t'+type_of_input] - dataset_neurons['threshold_t'+type_of_input])
    dataset_neurons['ΔV_thrp'] = abs (dataset_neurons['threshold_v' + type_of_input] - dataset_neurons['peak_v'+type_of_input])
    dataset_neurons['UpDown_ratio'] =  dataset_neurons['upstroke_downstroke_ratio' + type_of_input]


    features = list(dataset_neurons.columns[len(dataset_neurons.columns) - new_feat:len(dataset_neurons.columns)])

    features_ie = ['f_i_curve_slope', 'adaptation', 'Width', 'UpDown_width']

    # remove evident outliers
    dataset_neurons = dataset_neurons[dataset_neurons['upstroke_downstroke_ratio' + type_of_input] < 6]
    dataset_neurons = dataset_neurons[dataset_neurons['upstroke_downstroke_ratio' + type_of_input] >0]
    dataset_neurons = dataset_neurons[dataset_neurons['Width'] <0.01]
    dataset_neurons = dataset_neurons[dataset_neurons['Down_width'] < 0.01]
    dataset_neurons = dataset_neurons[dataset_neurons['Down_width'] >0]
    dataset_neurons = dataset_neurons[dataset_neurons['UpDown_width'] < 0.01]
    dataset_neurons = dataset_neurons[dataset_neurons['UpDown_width'] > 0]


    # let's find and remove the NaN for selected features
    feat_with_Nan=[]
    if dataset_type=='PatchSeqDataset':
        feat_with_Nan=['T-type Label']
    for feat in features:
        nan_count = dataset_neurons[feat].isna().sum()
        print("Nan for " + feat + ": " + str(nan_count))
        if nan_count!=0:
            feat_with_Nan.append(feat)

    if not len(feat_with_Nan):
        dataset_neurons_clean=dataset_neurons
    else:
        for feat in feat_with_Nan:
            dataset_neurons_clean = dataset_neurons[dataset_neurons[feat].notna()]

    feat_with_Nan=[]
    for feat in features_ie:
        nan_count = dataset_neurons_clean[feat].isna().sum()
        print("Nan for " + feat + ": " + str(nan_count))
        if nan_count!=0:
            feat_with_Nan.append(feat)

    if not len(feat_with_Nan):
        dataset_neurons_clean = dataset_neurons_clean
    else:
        for feat in feat_with_Nan:
            dataset_neurons_clean = dataset_neurons_clean[dataset_neurons_clean[feat].notna()]

    if dataset_type == 'PatchClampDataset':
        
        # save type of neurons for relevant rows
        type = dataset_neurons_clean['tag__dendrite_type']
        donor = dataset_neurons_clean['donor__species']
        region = dataset_neurons_clean['structure__id']
        specimen = dataset_neurons_clean['line_name']
        celltype_id=dataset_neurons_clean['specimen_id']
        transcript_id=[]

    if dataset_type == "PatchSeqDataset":
        
        #print("PatchSeqDataset pre-processing: ")
        
        #save type of neurons for relevant rows
        dataset_neurons_clean['T-type Label'].fillna('Unknown  ', inplace=True)
        specimen = dataset_neurons_clean['T-type Label'].map(lambda x: x.split(" ")[0])

        type = dataset_neurons_clean['dendrite_type']
        donor = dataset_neurons_clean['donor__species']
        #map each structure to an integer id for visualization
        region= (dataset_neurons_clean['structure'].astype('category')).cat.codes 
        transcript_id=dataset_neurons_clean['transcriptomics_sample_id']
        celltype_id=dataset_neurons_clean['cell_specimen_id']

    #print(list(dataset_neurons_clean.columns))
    print("Final number of cells after data cleaning "+str(len(dataset_neurons_clean.index)))

    # Histograms and kde
    
    # create a list of dataframe columns to use
    cols = ['UpDown_ratio', 'ΔV_deep', 'Slope_deep', 'ΔV_ratio']  
    sns.set(font_scale=2)


    for col in cols:
        plot=sns.displot(data=dataset_neurons_clean, x=col, kde=True, color='#824C71')
        plot.set_axis_labels(col, "Absolute Frequency")
        plt.savefig(os.path.join(path_for_images, 'histograms_feature_distributions_'+col+'.png'), dpi=300.0, bbox_inches='tight')
        # plt.show()

    #normalize features using z-score
    dataset_neurons_znorm=(dataset_neurons_clean[features]-dataset_neurons_clean[features].mean())/dataset_neurons_clean[features].std()

    return dataset_neurons_znorm, type, donor, region, specimen, transcript_id, celltype_id


def load_dataPatchSeqDataset(data_to_load='aspiny'):
    """
    Function to load PatchSeqDataset dataset from csv files.

    @param data_to_load: specify neurons to select ('full', 'aspiny', 'spiny', 'sparsely spiny') (string)

    @return dataset_filt: raw dataset considering only neurons from 'Mus musculus' selected according to 'data_to_load' (Pandas Dataframe)
    """
    # retrieve the data from the two csv files
    dataset=pd.read_csv(os.path.join("..", "data", "PatchSeqDataset", "PatchSeq_EP_features.csv"), delimiter=',')
    metadata=pd.read_csv(os.path.join("..", "data", "PatchSeqDataset", "PatchSeq_metadata.csv"), delimiter=',')

    # merge according to cell_specifimen_id removing duplicated columns
    cols_to_use = dataset.columns.difference(metadata.columns)
    dataset_ridotto= dataset[cols_to_use]
    dataset_ridotto.loc[:,('cell_specimen_id')]=dataset['cell_specimen_id']
    dfNew = pd.merge(dataset_ridotto, metadata, on='cell_specimen_id')

    # filter according to neuron type, if requested
    if data_to_load != 'full':
        dataset_filt = dfNew[dfNew['dendrite_type'] == data_to_load]
    else:
        dataset_filt = dfNew

    dataset_filt['donor__species']='Mus musculus'

    #if voltage values are in nanovolt, put them in microvolt for consistency with PatchClampDataset
    for feat in ['fast_trough_v_ramp', 'peak_v_ramp', 'slow_trough_v_ramp','threshold_v_ramp', 'trough_v_ramp']:
        dataset_filt.loc[abs(dataset_filt[feat])>100, feat]=dataset_filt.loc[abs(dataset_filt[feat])>100, feat]/1000000000

    return dataset_filt

def load_dataPatchClampDataset(data_to_load='aspiny'):
    """
      Function to load PatchClampDataset dataset from csv files.

      @param data_to_load: specify neurons to select ('full', 'aspiny', 'spiny', 'sparsely spiny') (string)

      @return dataset_filt: raw dataset considering only neurons from 'Mus musculus' selected according to 'data_to_load' (Pandas Dataframe)
    """
    #retrieve the data from the two csv files
    dataset_long=pd.read_csv(os.path.join("..", "data", "PatchClampDataset", "PatchClamp_metadata.csv"), delimiter=';')
    dataset_short=pd.read_csv(os.path.join("..", "data", "PatchClampDataset", "PatchClamp_EP_features.csv"), delimiter=',')

    #merge according to specifimen id removing duplicated columns
    cols_to_use = dataset_long.columns.difference(dataset_short.columns)
    dataset_long_reduced=dataset_long[cols_to_use]
    dataset_long_reduced.loc[:,('specimen_id')]=dataset_long['specimen_id']
    dfNew = pd.merge(dataset_short, dataset_long_reduced, on='specimen_id')

    # filter according to neuron type, if requested
    if data_to_load != 'full':
        dataset_filt = dfNew[dfNew['tag__dendrite_type'] == data_to_load]
    else:
        dataset_filt=dfNew

    #analyse only data from mice
    dataset_filt = dataset_filt[dataset_filt['donor__species'] == 'Mus musculus']

    return dataset_filt


def population_change_name(inputFile, outputFile):
    
    inFile = open(inputFile, "r")
    outFile = open(outputFile, "wt")

    pvalb_family = ["Pvalb-IRES-Cre", "Pvalb-T2A-Dre", "Nkx2-1-CreERT2", "Htr3a-Cre_NO152|Pvalb-T2A-Dre",
                    "Pvalb-T2A-FlpO|Vipr2-IRES2-Cre", "Chrna2-Cre_OE25|Pvalb-T2A-Dre", "Slc32a1-T2A-FlpO",
                    "Pvalb-T2A-CreERT2", "Slc32a1-T2A-FlpO|Vipr2-IRES2-Cre"]
    sst_family = ["Sst-IRES-Cre", "Sst-IRES-FlpO", "Chrna2-Cre_OE25", "Nos1-CreERT2|Sst-IRES-FlpO",
                  "Sst-IRES-Cre-197628.06.02.01"]
    htr3a_family = ["Ndnf-IRES2-dgCre", "Htr3a-Cre_NO152", "Chat-IRES-Cre", "Vip-IRES-Cre", "Chat-IRES-Cre-neo"]
    unsure_family = ["Oxtr-2A-Cre", "Oxtr-T2A-Cre", "Vipr2-IRES2-Cre", "Nos1-CreERT2", "Gad2-IRES-Cre", "Cux2-CreERT2",
                     "Rorb-IRES2-Cre", "Scnn1a-Tg2-Cre", "Gng7-Cre_KH71", "Nr5a1-Cre", "Scnn1a-Tg3-Cre", "Ntsr1-Cre",
                     "Rbp4-Cre_KL100", "Ctgf-T2A-dgCre", "Esr2-IRES2-Cre", "Ntsr1-Cre_GN220"]


    print("Generating files reporting Cre lines cell types prevalences...")
    for line in inFile:

        l_name = line.rstrip()
        l_name = l_name.split(",")[0]

        if l_name in pvalb_family:
            outFile.write(line.replace(l_name, "Pvalb"))

        elif l_name in sst_family:
            outFile.write(line.replace(l_name, "Sst"))

        elif l_name in htr3a_family:
            outFile.write(line.replace(l_name, "Vip/Lamp5"))

        elif l_name in unsure_family:
            outFile.write(line.replace(l_name, "Unsure"))

        else:
            outFile.write(line)

    inFile.close()
    outFile.close()


def cluster_population_lines(file):
    
    print("Analyzing cell type composition...")
    inf = open(file, "r")

    line = {}

    for riga in inf:
        l_name = riga.rstrip()
        l_name = l_name.split(",")[0]
        if l_name in line:
            line[l_name] = line[l_name] + 1
        else:
            line[l_name] = 1

    #print(line)
    print("There are ", len(line), " specimen lines")

    # TODO: esplicitare da dove abbiamo estratto queste liste
    pvalb_family = ["Pvalb-IRES-Cre", "Pvalb-T2A-Dre", "Nkx2-1-CreERT2", "Htr3a-Cre_NO152|Pvalb-T2A-Dre",
                    "Pvalb-T2A-FlpO|Vipr2-IRES2-Cre", "Chrna2-Cre_OE25|Pvalb-T2A-Dre", "Slc32a1-T2A-FlpO",
                    "Pvalb-T2A-CreERT2", "Slc32a1-T2A-FlpO|Vipr2-IRES2-Cre"]
    sst_family = ["Sst-IRES-Cre", "Sst-IRES-FlpO", "Chrna2-Cre_OE25", "Nos1-CreERT2|Sst-IRES-FlpO",
                  "Sst-IRES-Cre-197628.06.02.01"]
    htr3a_family = ["Ndnf-IRES2-dgCre", "Htr3a-Cre_NO152", "Chat-IRES-Cre", "Vip-IRES-Cre", "Chat-IRES-Cre-neo"]
    unsure_family = ["Oxtr-2A-Cre", "Oxtr-T2A-Cre", "Vipr2-IRES2-Cre", "Nos1-CreERT2", "Gad2-IRES-Cre", "Cux2-CreERT2",
                     "Rorb-IRES2-Cre", "Scnn1a-Tg2-Cre", "Gng7-Cre_KH71", "Nr5a1-Cre", "Scnn1a-Tg3-Cre", "Ntsr1-Cre",
                     "Rbp4-Cre_KL100", "Ctgf-T2A-dgCre", "Esr2-IRES2-Cre", "Ntsr1-Cre_GN220"]

    types = {"pvalb_family": [pvalb_family, 0], "sst_family": [sst_family, 0],
             "htr3a_family": [htr3a_family, 0], "unsure_family": [unsure_family, 0]}

    unnamed_count = 0
    for i in line:

        if i in pvalb_family:
            types["pvalb_family"][1] = types["pvalb_family"][1] + line[i]

        elif i in sst_family:
            types["sst_family"][1] = types["sst_family"][1] + line[i]

        elif i in htr3a_family:
            types["htr3a_family"][1] = types["htr3a_family"][1] + line[i]

        elif i in unsure_family:
            types["unsure_family"][1] = types["unsure_family"][1] + line[i]

        else:
            unnamed_count = unnamed_count + line[i]

    tot = unnamed_count
    for i in types:
        tot = tot + types[i][1]

    #for i in types:
        #print("# of cell types of ", i, " = ", types[i][1])
        #print("percentage : ", types[i][1] / tot * 100)

    #print("# of cell types of unnamed = ", unnamed_count)
    #print("percentage : ", unnamed_count / tot * 100)

    
    list_pop = [types["pvalb_family"][1], types["sst_family"][1], types["htr3a_family"][1], tot]
    inf.close()

    return list_pop


if __name__=='__main__':
    if (len(sys.argv)!=3):
        print("Wrong number of arguments! Format is: EP_analysis.py [dataset_type] [number_of_clusters] ")
    dataset_type=sys.argv[1]
    if sys.argv[2].isdigit():
        number_of_clusters=int(sys.argv[2])
    else:
        #will set number_of_clusters = max(silhouette), but here first set it to -1 as flag
        number_of_clusters=-1
    main(dataset_type, 'aspiny', 'ramp', number_of_clusters, 'CLUSTER')
