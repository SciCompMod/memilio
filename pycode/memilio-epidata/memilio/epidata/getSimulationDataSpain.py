import pandas as pd 
import os
import io
import requests

def fetch_population_data():
    download_url = 'https://servicios.ine.es/wstempus/js/es/DATOS_TABLA/67988?tip=AM&'
    req = requests.get(download_url)
    req.encoding = 'ISO-8859-1' 

    df = pd.read_json(io.StringIO(req.text))
    df = df[['MetaData', 'Data']]
    df['ID_Provincia'] = df['MetaData'].apply(lambda x: x[0]['Id'])
    df['Population'] = df['Data'].apply(lambda x: x[0]['Valor'])
    return df[['ID_Provincia', 'Population']]

def remove_islands(df):
    df = df[~df['ID_Provincia'].isin([51, 52, 8, 35, 38])]
    return df

def get_population_data():
    df = fetch_population_data()
    df = remove_islands(df)

    return df

def fetch_icu_data():
    download_url = 'https://www.sanidad.gob.es/areas/alertasEmergenciasSanitarias/alertasActuales/nCov/documentos/Datos_Capacidad_Asistencial_Historico_14072023.csv'
    req = requests.get(download_url)
    req.encoding = 'ISO-8859-1' 

    df = pd.read_csv(io.StringIO(req.text), sep=';')
    return df

def preprocess_icu_data(df):
    df_without_ventilator = df[df["Unidad"] == "U. Críticas SIN respirador"]
    df_with_ventilator = df[df["Unidad"] == "U. Críticas CON respirador"]

    df_icu = df_without_ventilator[['Fecha', 'ID_Provincia', 'Provincia', 'OCUPADAS_COVID19']].rename(columns={'OCUPADAS_COVID19': 'ICU'})
    df_icu_vent = df_with_ventilator[['Fecha', 'ID_Provincia', 'Provincia', 'OCUPADAS_COVID19']].rename(columns={'OCUPADAS_COVID19': 'ICU_ventilated'})

    df_merged = pd.merge(df_icu, df_icu_vent, on=['Fecha', 'ID_Provincia', 'Provincia'], how='outer')
    df_merged['Fecha'] = pd.to_datetime(df_merged['Fecha'], format='%d/%m/%Y')
    return df_merged

def get_icu_data():
    df = fetch_icu_data()
    df.rename(columns={'Cod_Provincia': 'ID_Provincia'}, inplace=True)
    df = remove_islands(df)
    df = preprocess_icu_data(df)

    return df

if __name__ == "__main__":

    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../data/Spain")

    df = get_population_data()
    df.to_json(os.path.join(data_dir, 'pydata/provincias_current_population.json'), orient='records')

    df = get_icu_data()
    df.to_json(os.path.join(data_dir, 'pydata/provincia_icu.json'), orient='records')