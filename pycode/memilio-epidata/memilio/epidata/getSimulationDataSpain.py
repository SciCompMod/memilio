import pandas as pd
import os
import io
import requests

from pyspainmobility import Mobility, Zones
import defaultDict as dd

from memilio.epidata import getDataIntoPandasDataFrame as gd


def fetch_population_data():
    download_url = 'https://servicios.ine.es/wstempus/js/es/DATOS_TABLA/67988?tip=AM&'
    req = requests.get(download_url)
    req.encoding = 'ISO-8859-1'

    df = pd.read_json(io.StringIO(req.text))
    df = df[['MetaData', 'Data']]
    df = df[df['MetaData'].apply(
        lambda x: x[0]['T3_Variable'] == 'Provincias')]
    df = df[df['MetaData'].apply(
        lambda x: x[1]['Nombre'] == 'Total')]
    df['ID_Provincia'] = df['MetaData'].apply(lambda x: x[0]['Id'])
    df['Population'] = df['Data'].apply(lambda x: x[0]['Valor'])
    return df[['ID_Provincia', 'Population']]


def remove_islands(df, column_labels=['ID_Provincia']):
    for label in column_labels:
        df = df[~df[label].isin([51, 52, 8, 35, 38])]
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
    df_icu = df[df["Unidad"] == "U. Críticas SIN respirador"]
    df_icu_vent = df[df["Unidad"] == "U. Críticas CON respirador"]

    df_icu = df_icu[['Fecha', 'ID_Provincia', 'OCUPADAS_COVID19']].rename(
        columns={'OCUPADAS_COVID19': 'ICU'})
    df_icu_vent = df_icu_vent[['Fecha', 'ID_Provincia', 'OCUPADAS_COVID19']].rename(
        columns={'OCUPADAS_COVID19': 'ICU_ventilated'})

    df_merged = pd.merge(df_icu, df_icu_vent, on=[
                         'Fecha', 'ID_Provincia'], how='outer')
    df_merged['Fecha'] = pd.to_datetime(
        df_merged['Fecha'], format='%d/%m/%Y').dt.strftime('%Y-%m-%d')
    df_merged.rename(columns={'Fecha': 'Date'}, inplace=True)

    return df_merged


def get_icu_data():
    df = fetch_icu_data()
    df.rename(columns={'Cod_Provincia': 'ID_Provincia'}, inplace=True)
    df = remove_islands(df)
    df = preprocess_icu_data(df)

    return df


def fetch_case_data():
    download_url = 'https://cnecovid.isciii.es/covid19/resources/casos_diagnostico_provincia.csv'
    req = requests.get(download_url)

    df = pd.read_csv(io.StringIO(req.text), sep=',')

    return df


def preprocess_case_data(df):
    df['provincia_iso'] = df['provincia_iso'].map(dd.Provincia_ISO_to_ID)
    df = df.rename(
        columns={'provincia_iso': 'ID_County', 'fecha': 'Date',
                 'num_casos': 'Confirmed'})[
        ['ID_County', 'Date', 'Confirmed']]
    # ensure correct types
    df['ID_County'] = pd.to_numeric(
        df['ID_County'], errors='coerce').astype('Int64')
    df = df[df['ID_County'].notna()].copy()
    df['ID_County'] = df['ID_County'].astype(int)
    return df


def get_case_data():
    df_raw = fetch_case_data()
    df = preprocess_case_data(df_raw)
    df = df.groupby(['ID_County', 'Date'], as_index=False)['Confirmed'].sum()
    df = df.sort_values(['ID_County', 'Date'])
    df['Confirmed'] = df.groupby('ID_County')['Confirmed'].cumsum()
    return df


def download_mobility_data(data_dir, start_date, end_date, level):
    mobility_data = Mobility(version=2, zones=level,
                             start_date=start_date, end_date=end_date, output_directory=data_dir)
    mobility_data.get_od_data()


def preprocess_mobility_data(df, data_dir):
    df.drop(columns=['hour', 'trips_total_length_km'], inplace=True)
    if not os.path.exists(os.path.join(data_dir, 'poblacion.csv')):
        download_url = 'https://movilidad-opendata.mitma.es/zonificacion/poblacion.csv'
        req = requests.get(download_url)
        with open(os.path.join(data_dir, 'poblacion.csv'), 'wb') as f:
            f.write(req.content)

    zonification_df = pd.read_csv(os.path.join(data_dir, 'poblacion.csv'), sep='|')[
        ['municipio', 'provincia', 'poblacion']]
    poblacion = zonification_df.groupby('provincia')[
        'poblacion'].sum()
    zonification_df.drop_duplicates(
        subset=['municipio', 'provincia'], inplace=True)
    municipio_to_provincia = dict(
        zip(zonification_df['municipio'], zonification_df['provincia']))
    df['id_origin'] = df['id_origin'].map(municipio_to_provincia)
    df['id_destination'] = df['id_destination'].map(municipio_to_provincia)

    df.query('id_origin != id_destination', inplace=True)
    df = df.groupby(['date', 'id_origin', 'id_destination'],
                    as_index=False)['n_trips'].sum()

    df['n_trips'] = df.apply(
        lambda row: row['n_trips'] / poblacion[row['id_origin']], axis=1)

    df = df.groupby(['id_origin', 'id_destination'],
                    as_index=False)['n_trips'].mean()

    df = remove_islands(df, ['id_origin', 'id_destination'])

    return df


def get_mobility_data(data_dir, start_date='2022-08-01', end_date='2022-08-31', level='municipios'):
    filename = f'Viajes_{level}_{start_date}_{end_date}_v2.parquet'
    if not os.path.exists(os.path.join(data_dir, filename)):
        print(
            f"File {os.path.exists(os.path.join(data_dir, filename))} does not exist. Downloading mobility data...")
        download_mobility_data(data_dir, start_date, end_date, level)

    df = pd.read_parquet(os.path.join(data_dir, filename))
    df = preprocess_mobility_data(df, data_dir)

    return df


if __name__ == "__main__":

    data_dir = os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "../../../../data/Spain")
    pydata_dir = os.path.join(data_dir, 'pydata')
    os.makedirs(pydata_dir, exist_ok=True)
    mobility_dir = os.path.join(data_dir, 'mobility')
    os.makedirs(mobility_dir, exist_ok=True)

    df = get_population_data()
    df.to_json(os.path.join(
        pydata_dir, 'provincias_current_population.json'), orient='records')

    df = get_icu_data()
    # rename to match expected column name in parameters_io.h
    if 'ID_Provincia' in df.columns:
        df.rename(columns={'ID_Provincia': 'ID_County'}, inplace=True)
    # and should also be int
    if 'ID_County' in df.columns:
        df['ID_County'] = pd.to_numeric(df['ID_County'], errors='coerce')
        df = df[df['ID_County'].notna()].copy()
        df['ID_County'] = df['ID_County'].astype(int)
    df.to_json(os.path.join(pydata_dir, 'provincia_icu.json'), orient='records')

    df = get_case_data()
    # same for case data
    if 'ID_Provincia' in df.columns:
        df.rename(columns={'ID_Provincia': 'ID_County'}, inplace=True)
    if 'ID_County' in df.columns:
        df['ID_County'] = pd.to_numeric(df['ID_County'], errors='coerce')
        df = df[df['ID_County'].notna()].copy()
        df['ID_County'] = df['ID_County'].astype(int)
    df.to_json(os.path.join(
        pydata_dir, 'cases_all_pronvincias.json'), orient='records')

    df = get_mobility_data(mobility_dir)
    matrix = df.pivot(index='id_origin',
                      columns='id_destination', values='n_trips').fillna(0)

    gd.write_dataframe(matrix, data_dir, 'commuter_mobility', 'txt', {
                       'sep': ' ', 'index': False, 'header': False})
