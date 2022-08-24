# Import packages
import os
import pickle
import pandas as pd
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt
# Ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
#%%
# DOWNLOAD DEPMAP DATA


class DataDownload:

    """
    Class for downloading DepMap data from the specified URL, wrangling, cleaning and saving as a new file
    """

    def __init__(self,
                 col_id='DepMap_ID',
                 cellname='stripped_cell_line_name',
                 rename='CCLE_Name',
                 foldername="foldername",
                 filename="csvfilename"
                 ):
        """
        Constructor for DataDownload Class

        :param col_id: (STR) Name of unique identifier column
        :param cellname: (STR) Name of stripped cell line name column
        :param rename: (STR) Preferred column name to rename stripped cell line
        :param foldername: (STR) Name of directory to create to store data file
        :param filename: (STR) Name of data file
        """

        self.col_id = col_id
        self.cellname = cellname
        self.rename = rename
        self.foldername = foldername
        self.filename = filename

    @staticmethod
    def data_download(url):
        """
        Download data file from specified URL and return as a pandas DF

        :param url: (STR) URL of data download
        :return df: (DF) Pandas DF of downloaded data
        """
        df = pd.read_csv(url)
        return df

    def data_organise(self, df_sampleinfo, df_geneeffect):
        """
        From df_sampleinfo, keep columns: DepMap_ID and stripped_cell_line_name
        Rename stripped_cell_line_name as CCLE_Name.
        Pandas.merge reduced dfs on DepMap_ID column of df_geneeffect

        :param df_sampleinfo: (DF) Pandas DF of sample info DepMap datafile
        :param df_geneeffect: (DF) Pandas DF of gene_effect DepMap datafile
        :return df_reindex: (DF) Merged DF of df_sampleinfo and df_geneeffect
        """

        # Drop columns and merge based on df_geneeffect DepMap ID
        df_sampleinfo = df_sampleinfo[[self.col_id,
                                       self.cellname
                                       ]]
        df_sampleinfo = df_sampleinfo.rename(columns={self.cellname: self.rename})
        df_merged = pd.merge(df_geneeffect,
                             df_sampleinfo,
                             on=self.col_id,
                             how='left'
                             )
        # Reset index of df_merged
        df_reindex = df_merged.drop([self.col_id],
                                    axis=1
                                    )
        df_reindex = df_reindex.set_index([self.rename])
        df_reindex.columns = df_reindex.columns.str.split(' ', 1).str.get(0)
        return df_reindex

    @staticmethod
    def fill_nan(df):
        """
        Identify missing data points and impute values using SKLearn SimpleImputer.
        Impute the median from the column of the NaN.
        Median imputation due to non-normally distributed skewed data.

        :param df: (DF) Wrangled and merged DF of sample_info and gene_effect data
        :return df_final: (DF) DF with NaN values imputed using median value
        """

        nan_sum = df.isna().sum().sum()
        imputer = SimpleImputer(strategy='median')
        # Only impute values if number of missing values > 0
        if nan_sum > 0:
            imputer.fit(df)
            median_value = imputer.transform(df)
            df_final = pd.DataFrame(median_value,
                                    columns=df.columns,
                                    index=df.index
                                    )
        else:
            df_final = df
        print(f"{nan_sum} missing values imputed with median from each column")
        return df_final

    @staticmethod
    def load_data(file, index, directory):
        """
        Load CSV file as a pandas DF

        :param file: (STR) Name of file
        :param index: (INT) Column of CSV file to use as DF index
        :param directory: (STR) Name of directory where CSV file is saved
        :return: (DF) Pandas DF of CSV file
        """

        path = os.path.join(directory,
                            file
                            )
        return pd.read_csv(path, index_col=index)

    def save_df(self, df):
        """
        Save DF  as a CSV file in project directory

        :param df: (DF) Pandas DF to save as CSV file
        :return: (CSV) DF saved as CSV
        """

        project_root_dir = ".."
        files = os.path.join(project_root_dir,
                             self.foldername
                             )
        os.makedirs(files,
                    exist_ok=True
                    )
        df.to_csv(os.path.join(files,
                               self.filename + '.csv'
                               )
                  )

    def clean_and_save(self, url_sampleinfo, url_geneeffect):
        """
        Calling function to wrangle, merge and clean final DepMap data frame and save as CSV file

        :param url_sampleinfo: (STR) URL of sample_info data file
        :param url_geneeffect: (STR) URL of gene_effect data file
        :return: (DF) Final DF of DepMap cell-fitness data saved as CSV file
        """
        # Download data as pandas df
        print('Downloading data from URLs')
        df_sampleinfo = self.data_download(url_sampleinfo)
        df_geneeffect = self.data_download(url_geneeffect)
        print('-'*20)
        # Organise and merge dataframes
        print('Organising and merging dataframes')
        df_organised = self.data_organise(df_sampleinfo, df_geneeffect)
        print('-'*20)
        # Fill NaN with imputed values
        print('Filling NaN values')
        df_filled = self.fill_nan(df_organised)
        print('-'*20)
        # Save df
        print('Saving clean data file')
        self.save_df(df_filled)
        print('-'*20)
        print('Complete')
        print('-'*20)

#%%
# SAVE RESULTS


class SaveLoad:
    """
    Class for saving and loading pickle files and figures
    """

    def __init__(self,
                 picklefolder="pickle_files",
                 figfolder="Figures"):
        """
        Constructor for SaveLoad Class

        :param picklefolder: (STR) Name of folder to save pickle file
        :param figfolder: (STR) Name of folder to save figures
        """

        self.picklefolder = picklefolder
        self.figfolder = figfolder

    def save_dict_pickle(self, dictex, filename):
        """
        Create a new directory to then save a dictionary as a pickle file

        :param dictex: (DICT) Dictionary to export as a pickle file
        :param filename: (STR) Name of file, should not include '.pickle'
        :return: (PICKLE) Dictionary saved as pickle file
        """

        # Create a folder to save the pickle files
        project_root_dir = "."
        directory = os.path.join(project_root_dir,
                                 self.picklefolder
                                 )
        os.makedirs(directory,
                    exist_ok=True
                    )

        # Open directory to save dictionary as pickle
        pickle_out = open(os.path.join(directory, f'{filename}.pickle'), 'wb')
        # Save dictionary as pickle file
        pickle.dump(dictex, pickle_out)
        pickle_out.close()

    @staticmethod
    def load_dict_pickle(file):
        """
        Load a pickle file as a python ductionary

        :param file: (STR) Path of file to load to include '.pickle'
        :return new_dict: (DICT) Pickle file in dictionary format
        """

        pickle_in = open(f'{file}', 'rb')
        new_dict = pickle.load(pickle_in)
        return new_dict

    def save_fig(self, fig_id, tight_layout=False, fig_extension="pdf", resolution=300):
        """
        Save a figure in zArchive2 directory. Default format is PDF

        :param fig_id: (STR) Figure name
        :param tight_layout: (BOOLEAN) If true save figure as tight layout, default=False
        :param fig_extension: (STR) Figure file format, default=PDF
        :param resolution: (INT) Figure resolution, default=300
        """

        # Create a folder to save the figure
        project_root_dir = ".."
        directory = os.path.join(project_root_dir,
                                 self.figfolder
                                 )
        os.makedirs(directory,
                    exist_ok=True
                    )
        print("Saving figure", fig_id)
        # Save the figure, if tight layout true save as tight layout
        if tight_layout:
            plt.tight_layout()
        plt.savefig(os.path.join(directory, fig_id + "." + fig_extension), format=fig_extension, dpi=resolution)
