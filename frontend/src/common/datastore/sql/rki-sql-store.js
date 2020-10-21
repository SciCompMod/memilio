import axios from 'axios';
import SQLDatastore, {SQLColumn} from './sql-datastore';
import SQLType from './SQLType';
import {stateIdFromCountyId} from '../../utils';

// TODO: Use multiple tables and JOINs?
/**
 * This class manages all RKI related data.
 */
class RKISQLStore extends SQLDatastore {
  /** @enum */
  Level = {
    GERMANY: 0,
    STATE: 1,
    COUNTY: 2,
  };

  /** @enum */
  Type = {
    ALL: 0,
    GENDERS: 1,
    AGE_GROUP: 2,
  };

  constructor() {
    super();
    this.table = {
      name: 'rki',
      columns: [
        new SQLColumn('level', SQLType.INT2, true),

        new SQLColumn('stateId', SQLType.INT),
        new SQLColumn('countyId', SQLType.INT),

        new SQLColumn('name', SQLType.TEXT),

        new SQLColumn('date', SQLType.BIGINT, true),

        new SQLColumn('confirmed', SQLType.INT, true),
        new SQLColumn('recovered', SQLType.INT, true),
        new SQLColumn('deaths', SQLType.INT, true),
      ],
    };

    this.create(this.table);
  }

  async populate() {
    let response = await axios.get('assets/all_state_rki.json');
    let data = response.data;
    let dataTransformed = [];

    for (const datum of data) {
      dataTransformed.push([
        this.Level.STATE,
        datum.ID_State,
        null,
        datum.State,
        datum.Date,
        datum.Confirmed,
        datum.Recovered,
        datum.Deaths,
      ]);
    }

    this.putAll(this.table.name, dataTransformed);

    response = await axios.get('assets/all_county_rki.json');

    data = response.data;
    dataTransformed = [];

    for (const datum of data) {
      dataTransformed.push([
        this.Level.COUNTY,
        stateIdFromCountyId(datum.ID_County),
        datum.ID_County,
        datum.County,
        datum.Date,
        datum.Confirmed,
        datum.Recovered,
        datum.Deaths,
      ]);
    }

    // We can't add 25000 entries at the same time due to memory limitations, so we slice them up into smaller batches.
    const batchSize = 2000;
    for (let i = 0; i < dataTransformed.length; i += batchSize) {
      const subArray = dataTransformed.slice(i, i + batchSize);
      await this.putAll(this.table.name, subArray);
    }
  }

  /**
   * Returns all entries for Germany ordered by date.
   * @return {Promise<Array<{date: number, confirmed: number, recovered: number, deaths: number}>>}
   */
  async getAllGermany() {
    // const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.GERMANY} ORDER BY date`;
    const query = `SELECT 
      date, 
      SUM(confirmed) as confirmed, 
      SUM(recovered) as recovered, 
      SUM(deaths) as deaths 
    FROM ${this.table.name}
    WHERE level=${this.Level.STATE}
    GROUP BY date
    ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  /**
   * Returns all entries for the given state.
   * @param stateId{number}
   * @param columns{string} Comma separated list of column names.
   * @return {Promise<Array<{date: number, confirmed: number, recovered: number, deaths: number}>>}
   */
  async getAllState(stateId, columns = 'date, confirmed, recovered, deaths') {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.STATE} AND stateId='${stateId}' ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  /**
   * Returns all entries for the given state in the specified date range.
   * @param stateId{number}
   * @param range{{start: number, end: number}}
   * @param columns{string} Comma separated list of column names.
   * @return {Promise<Array<{date: number, confirmed: number, recovered: number, deaths: number}>>}
   */
  async getAllStateInRange(stateId, range, columns = 'date, confirmed, recovered, deaths') {
    const query = `SELECT 
      ${columns} 
    FROM ${this.table.name} 
    WHERE level=${this.Level.STATE} 
      AND stateId='${stateId}' 
      AND (date BETWEEN ${range.start} AND ${range.end}) 
    ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  /**
   * Returns all entries for the given county.
   * @param countyId{number}
   * @param columns{string} Comma separated list of column names.
   * @return {Promise<Array<{date: number, confirmed: number, recovered: number, deaths: number}>>}
   */
  async getAllCounty(countyId, columns = 'date, confirmed, recovered, deaths') {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.COUNTY} AND countyId='${countyId}' ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  /**
   * Gets all entries for all states.
   * @param columns{string} Comma separated list of column names.
   * @return {Promise<Array<{date: number, stateId: number, confirmed: number, recovered: number, deaths: number}>>}
   */
  async getAllStates(columns = 'date, stateId, confirmed, recovered, deaths') {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.STATE} ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  /**
   * Gets all entries in the specified date range for all states.
   * @param range{{start: number, end: number}}
   * @param columns{string} Comma separated list of column names.
   * @return {Promise<Array<{date: number, stateId: number, confirmed: number, recovered: number, deaths: number}>>}
   */
  async getAllStatesInRange(range, columns = 'date, stateId, confirmed, recovered, deaths') {
    const query = `SELECT 
      ${columns} 
    FROM ${this.table.name} 
    WHERE level=${this.Level.STATE}  
      AND (date BETWEEN ${range.start} AND ${range.end}) 
    ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  /**
   * Returns all entries for counties in the given state.
   * @param stateId{number}
   * @param columns{string} Comma separated list of column names.
   * @return {Promise<Array<{date: number, countyId: number, confirmed: number, recovered: number, deaths: number}>>}
   */
  async getAllCountiesOfState(stateId, columns = 'date, countyId, confirmed, recovered, deaths') {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.COUNTY} AND stateId='${stateId}' ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  /**
   * Returns entries in the given date range for all counties.
   * @param range{{start: number, end: number}}
   * @param columns{string} Comma separated list of column names.
   * @return {Promise<Array<{date: number, countyId: number, confirmed: number, recovered: number, deaths: number}>>}
   */
  async getAllCountiesInRange(range, columns = 'date, countyId, confirmed, recovered, deaths') {
    const query = `SELECT 
      ${columns} 
    FROM ${this.table.name} 
    WHERE level=${this.Level.COUNTY}  
      AND (date BETWEEN ${range.start} AND ${range.end}) 
    ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  /**
   * Returns the minimum and maximum time of existing data.
   * @return {Promise<{start: number, end: number}>}
   */
  async getTimeBounds() {
    const query = `SELECT MIN(date) as start, MAX(date) as end FROM ${this.table.name} WHERE level=${this.Level.STATE}`;
    return this.resultToArray(await this.execQuery(query))[0];
  }
}

const RKIStore = new RKISQLStore();
export default RKIStore;
