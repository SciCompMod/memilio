import pandas as pd
import os
import io
import requests
import numpy as np

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
    df['ID_Provincia'] = df['MetaData'].apply(
        lambda x: dd.provincia_id_map_census[x[0]['Id']])
    df['Population'] = df['Data'].apply(lambda x: x[0]['Valor'])
    return df[['ID_Provincia', 'Population']]


def remove_islands(df, column_labels=['ID_Provincia']):
    for label in column_labels:
        df = df[~df[label].isin([530, 630, 640, 701, 702])]
    return df


def get_population_data():
    df = fetch_population_data()
    df = df.sort_values(by=['ID_Provincia'])

    age_groups = ["<3 years", "3-5 years", "6-14 years", "15-17 years",
                  "18-24 years", "25-29 years", "30-39 years", "40-49 years",
                  "50-64 years", "65-74 years", ">74 years"]

    # Move Population into the first age group, ensure integer type
    df[age_groups[0]] = pd.to_numeric(
        df['Population'], errors='coerce').fillna(0).astype(int)
    # Initialize the remaining age groups to zero
    for g in age_groups[1:]:
        df[g] = 0

    df = remove_islands(df)

    return df


def fetch_icu_data():
    # https://www.sanidad.gob.es/areas/alertasEmergenciasSanitarias/alertasActuales/nCov/capacidadAsistencial.htm
    download_url = 'https://www.sanidad.gob.es/areas/alertasEmergenciasSanitarias/alertasActuales/nCov/documentos/Datos_Capacidad_Asistencial_Historico_14072023.csv'
    req = requests.get(download_url)
    req.encoding = 'ISO-8859-1'

    df = pd.read_csv(io.StringIO(req.text), sep=';')
    return df


def preprocess_icu_data(df, moving_average):
    df_icu = df[df["Unidad"] == "U. Críticas SIN respirador"]
    df_icu_vent = df[df["Unidad"] == "U. Críticas CON respirador"]

    df_icu = df_icu[['Fecha', 'ID_Provincia', 'OCUPADAS_COVID19']].rename(
        columns={'OCUPADAS_COVID19': 'ICU'})
    df_icu_vent = df_icu_vent[['Fecha', 'ID_Provincia', 'OCUPADAS_COVID19']].rename(
        columns={'OCUPADAS_COVID19': 'ICU_ventilated'})

    df_merged = pd.merge(df_icu, df_icu_vent, on=[
                         'Fecha', 'ID_Provincia'], how='outer')
    # df_merged['ICU'] = df_merged['ICU'].fillna(
    #     0) + df_merged['ICU_ventilated'].fillna(0)
    df_merged['Fecha'] = pd.to_datetime(
        df_merged['Fecha'], format='%d/%m/%Y').dt.strftime('%Y-%m-%d')
    df_merged.rename(columns={'Fecha': 'Date'}, inplace=True)

    df_merged = df_merged.sort_values(by=['ID_Provincia', 'Date'])

    if moving_average > 0:
        # Calculate 7-day moving average for ICU and ICU_ventilated
        df_merged['ICU'] = df_merged.groupby('ID_Provincia')['ICU'].transform(
            lambda x: x.rolling(window=moving_average, min_periods=1).mean())
        df_merged['ICU_ventilated'] = df_merged.groupby('ID_Provincia')['ICU_ventilated'].transform(
            lambda x: x.rolling(window=moving_average, min_periods=1).mean())
    # print provinces where ICU or ICU_ventilated are zero for all dates
    zero_icu = df_merged.groupby('ID_Provincia')['ICU'].apply(
        lambda x: x.fillna(0).eq(0).all())
    zero_icu_vent = df_merged.groupby('ID_Provincia')['ICU_ventilated'].apply(
        lambda x: x.fillna(0).eq(0).all())
    ids = sorted(set(zero_icu[zero_icu].index.tolist(
    ) + zero_icu_vent[zero_icu_vent].index.tolist()))
    if ids:
        print("ID_Provincia with ICU or ICU_ventilated zero for all dates:", ids)
        for pid in ids:
            icu_zero = bool(zero_icu.get(pid, False))
            icu_vent_zero = bool(zero_icu_vent.get(pid, False))
            print(
                f"ID_Provincia {pid}: ICU zero={icu_zero}, ICU_ventilated zero={icu_vent_zero}")

    return df_merged


def get_icu_data(moving_average=0):
    df = fetch_icu_data()
    df.rename(columns={'Cod_Provincia': 'ID_Provincia'}, inplace=True)
    # ensure numeric, drop non-numeric and increment all IDs by 1
    df['ID_Provincia'] = df['ID_Provincia'].map(dd.provincia_id_map)
    df = remove_islands(df)
    df = preprocess_icu_data(df, moving_average)

    return df


def fetch_case_data():
    # https://datos.gob.es/es/catalogo/e05070101-evolucion-de-enfermedad-por-el-coronavirus-covid-19
    download_url = 'https://cnecovid.isciii.es/covid19/resources/casos_diagnostico_provincia.csv'
    req = requests.get(download_url)

    df = pd.read_csv(io.StringIO(req.text), sep=',',
                     keep_default_na=False, na_values=[])

    return df


def preprocess_case_data(df, moving_average):
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

    df = df.sort_values(by=['ID_County', 'Date'])
    # df['Deaths'] = 0
    # df['Recovered'] = 0

    if moving_average > 0:
        # Calculate 7-day moving average
        df['ICU'] = df.groupby('ID_County')['Confirmed'].transform(
            lambda x: x.rolling(window=moving_average, min_periods=1).mean())

    return df


def get_case_data(moving_average=0):
    df_raw = fetch_case_data()
    df = preprocess_case_data(df_raw, moving_average)
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
    zonification_df['provincia'] = zonification_df['provincia'].map(
        dd.provincia_id_map)
    zonification_df.dropna(inplace=True, subset=['provincia'])
    poblacion = zonification_df.groupby('provincia')[
        'poblacion'].sum()
    zonification_df.drop_duplicates(
        subset=['municipio', 'provincia'], inplace=True)
    municipio_to_provincia = dict(
        zip(zonification_df['municipio'], zonification_df['provincia']))
    df['id_origin'] = df['id_origin'].map(
        municipio_to_provincia)
    df['id_destination'] = df['id_destination'].map(
        municipio_to_provincia)

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
            f"File {os.path.join(data_dir, filename)} does not exist. Downloading mobility data...")
        download_mobility_data(data_dir, start_date, end_date, level)

    df = pd.read_parquet(os.path.join(data_dir, filename))
    df = preprocess_mobility_data(df, data_dir)

    return df


def get_comunidadid_to_provinciaids_map():
    """ Creates a hash map from state IDs to lists of county IDs

    :param merge_eisenach: Default: True] Defines whether the counties
        'Wartburgkreis' and 'Eisenach' are listed separately or combined
        as one entity 'Wartburgkreis'.
    :param zfill: Default: False]. Defines whether or not all IDs are returned
        as zero-filled strings. By default, integer maps are returned.
    :returns: State IDs to lists of county IDs map

    """
    provincia_ids = [k for k in dd.Provincias.keys() if k not in {
        530, 630, 640, 701, 702}]
    comunidad_ids = list(sorted({k // 10 for k in provincia_ids}))
    comunidad_to_provincia_table = [[] for i in range(len(comunidad_ids))]

    for id in provincia_ids:
        comunidad_to_provincia_table[comunidad_ids.index(id // 10)].append(id)

    return dict(zip(comunidad_ids, comunidad_to_provincia_table))


def createComunidadesMobility(directory, mobility_file):
    mobility_matrix = pd.read_csv(
        os.path.join(directory + mobility_file + '.txt'),
        sep=' ', header=None)

    provinciaIDs = [k for k in dd.Provincias.keys() if k not in {
        530, 630, 640, 701, 702}]
    if (len(mobility_matrix.index) == len(provinciaIDs)) and (len(mobility_matrix.columns) == len(provinciaIDs)):
        # get county and state IDs
        comunidadIDs = list(sorted({k // 10 for k in provinciaIDs}))
        # get state ID to county ID map
        comunidadID_to_provinciaID = get_comunidadid_to_provinciaids_map()
        # initialize federal states mobility matrix
        mobility_matrix_comunidades = np.zeros(
            (len(comunidadID_to_provinciaID), len(comunidadID_to_provinciaID)))

        # iterate over state_to_county map and replace IDs by numbering 0, ..., n
        comunidad_indices = []
        provincia_indices = []
        for comunidad, privincias in comunidadID_to_provinciaID.items():
            comunidad_indices.append(comunidadIDs.index(comunidad))
            provincia_indices.append(
                np.array([provinciaIDs.index(provincia) for provincia in privincias]))

        comunidad_to_provincia = dict(
            zip(comunidad_indices, provincia_indices))
        # iterate over all comunidades
        for comunidad, provincias in comunidad_to_provincia.items():
            # iterate over all neighbors
            for comunidad_neighb, provincias_neighb in comunidad_to_provincia.items():

                if comunidad != comunidad_neighb:

                    mobility_matrix_comunidades[comunidad, comunidad_neighb] = mobility_matrix.iloc[provincias, provincias_neighb].sum(
                        axis=0).sum()
                    mobility_matrix_comunidades[comunidad_neighb, comunidad] = mobility_matrix.iloc[provincias_neighb, provincias].sum(
                        axis=1).sum()

        mobility_matrix_comunidades = pd.DataFrame(mobility_matrix_comunidades)
        gd.write_dataframe(
            mobility_matrix_comunidades, directory, mobility_file + '_comunidades', 'txt',
            param_dict={'sep': ' ', 'header': None, 'index': False})
        return mobility_matrix_comunidades
    else:
        return pd.DataFrame()


def aggregate_to_comunidades(df, column_labels=['ID_County']):
    df_aggregated = df.copy()
    for label in column_labels:
        df_aggregated[label] = df[label] // 10
    if 'Date' in df_aggregated.columns:
        df_aggregated = df_aggregated.groupby(
            column_labels + ['Date'], as_index=False).sum()
    else:
        df_aggregated = df_aggregated.groupby(
            column_labels, as_index=False).sum()
    return df_aggregated


if __name__ == "__main__":

    moving_average = 0

    data_dir = os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "../../../../data/Spain")
    pydata_dir = os.path.join(data_dir, 'pydata')
    os.makedirs(pydata_dir, exist_ok=True)
    mobility_dir = os.path.join(data_dir, 'mobility/')
    os.makedirs(mobility_dir, exist_ok=True)

    df = get_population_data()
    if 'ID_Provincia' in df.columns:
        df.rename(columns={'ID_Provincia': 'ID_County'}, inplace=True)
    df.to_json(os.path.join(
        pydata_dir, 'provincias_current_population.json'), orient='records')
    df_agg = aggregate_to_comunidades(df)
    df_agg.rename(columns={'ID_County': 'ID_State'}, inplace=True)
    df_agg.to_json(os.path.join(
        pydata_dir, 'comunidades_current_population.json'), orient='records')
    df_agg['ID_Country'] = 0
    df_agg = df_agg.drop(
        columns=['ID_State']).groupby(
        'ID_Country', as_index=True).sum()
    df_agg.to_json(os.path.join(
        pydata_dir, 'spain_current_population.json'), orient='records')

    df = get_icu_data(moving_average=moving_average)
    # rename to match expected column name in parameters_io.h
    if 'ID_Provincia' in df.columns:
        df.rename(columns={'ID_Provincia': 'ID_County'}, inplace=True)
    # and should also be int
    if 'ID_County' in df.columns:
        df['ID_County'] = pd.to_numeric(df['ID_County'], errors='coerce')
        df = df[df['ID_County'].notna()].copy()
        df['ID_County'] = df['ID_County'].astype(int)
    if moving_average > 0:
        df.to_json(os.path.join(pydata_dir, 'provincia_icu_ma' +
                   str(moving_average) + '.json'), orient='records')
        df_agg = aggregate_to_comunidades(df)
        df_agg.rename(columns={'ID_County': 'ID_State'}, inplace=True)
        df_agg.to_json(os.path.join(pydata_dir, 'comunidades_icu_ma' +
                                    str(moving_average) + '.json'), orient='records')
        df_agg['ID_Country'] = 0
        df_agg = df_agg.drop(
            columns=['ID_State']).groupby(
            ['ID_Country', 'Date'], as_index=False).sum()
        df_agg.drop(columns=['ID_Country'], inplace=True)
        df_agg.to_json(os.path.join(pydata_dir, 'spain_icu_ma' +
                                    str(moving_average) + '.json'), orient='records')
    else:
        df.to_json(os.path.join(pydata_dir, 'provincia_icu.json'),
                   orient='records')
        df_agg = aggregate_to_comunidades(df)
        df_agg.rename(columns={'ID_County': 'ID_State'}, inplace=True)
        df_agg.to_json(os.path.join(
            pydata_dir, 'comunidades_icu.json'), orient='records')
        df_agg['ID_Country'] = 0
        df_agg = df_agg.drop(
            columns=['ID_State']).groupby(
            ['ID_Country', 'Date'], as_index=False).sum()
        df_agg.drop(columns=['ID_Country'], inplace=True)
        df_agg.to_json(os.path.join(
            pydata_dir, 'spain_icu.json'), orient='records')

    df = get_case_data(moving_average=moving_average)
    # same for case data
    if 'ID_Provincia' in df.columns:
        df.rename(columns={'ID_Provincia': 'ID_County'}, inplace=True)
    if 'ID_County' in df.columns:
        df['ID_County'] = pd.to_numeric(df['ID_County'], errors='coerce')
        df = df[df['ID_County'].notna()].copy()
        df['ID_County'] = df['ID_County'].astype(int)
    if moving_average > 0:
        df.to_json(os.path.join(pydata_dir, 'cases_all_pronvincias_ma' +
                   str(moving_average) + '.json'), orient='records')
        df_agg = aggregate_to_comunidades(df, column_labels=['ID_County'])
        df_agg.rename(columns={'ID_County': 'ID_State'}, inplace=True)
        df_agg.to_json(os.path.join(pydata_dir, 'cases_all_comunidades_ma' +
                                    str(moving_average) + '.json'), orient='records')
        df_agg['ID_Country'] = 0
        df_agg = df_agg.drop(
            columns=['ID_State']).groupby(
            ['ID_Country', 'Date'], as_index=False).sum()
        df_agg.drop(columns=['ID_Country'], inplace=True)
        df_agg.to_json(os.path.join(pydata_dir, 'cases_all_spain_ma' +
                                    str(moving_average) + '.json'), orient='records')
    else:
        df.to_json(os.path.join(
            pydata_dir, 'cases_all_pronvincias.json'), orient='records')
        df_agg = aggregate_to_comunidades(df, column_labels=['ID_County'])
        df_agg.rename(columns={'ID_County': 'ID_State'}, inplace=True)
        df_agg.to_json(os.path.join(
            pydata_dir, 'cases_all_comunidades.json'), orient='records')
        df_agg['ID_Country'] = 0
        df_agg = df_agg.drop(
            columns=['ID_State']).groupby(
            ['ID_Country', 'Date'], as_index=False).sum()
        df_agg.drop(columns=['ID_Country'], inplace=True)
        df_agg.to_json(os.path.join(
            pydata_dir, 'cases_all_spain.json'), orient='records')

    df = get_mobility_data(mobility_dir)
    matrix = df.pivot(index='id_origin',
                      columns='id_destination', values='n_trips').fillna(0)

    gd.write_dataframe(matrix, mobility_dir, 'commuter_mobility', 'txt', {
                       'sep': ' ', 'index': False, 'header': False})

    createComunidadesMobility(mobility_dir, 'commuter_mobility')
