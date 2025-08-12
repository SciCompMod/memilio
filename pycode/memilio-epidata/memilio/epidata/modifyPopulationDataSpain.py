import pandas as pd 
import os

def read_population_data(file):
    df = pd.read_json(file)
    df = df[['MetaData', 'Data']]
    df['ID_Provincia'] = df['MetaData'].apply(lambda x: x[0]['Id'])
    df['Population'] = df['Data'].apply(lambda x: x[0]['Valor'])
    return df[['ID_Provincia', 'Population']]

def remove_islands(df):
    df = df[~df['ID_Provincia'].isin([51, 52, 8, 35, 38])]
    return df

if __name__ == "__main__":

    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../data/Spain")

    df = read_population_data(os.path.join(data_dir, 'pydata/67988.json'))
    df = remove_islands(df)
    df.to_json(os.path.join(data_dir, 'pydata/provincias_current_population.json'), orient='records', force_ascii=False)