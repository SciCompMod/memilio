import DB from './sql-database';
import axios from 'axios';
import SQLDatastore, {SQLColumn} from './sql-datastore';
import SQLType from './SQLType';
import {stateIdFromCountyId} from './utils';

// TODO: Use multiple tables and JOINs?
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
    super({
      name: 'rki',
      columns: [
        new SQLColumn('level', SQLType.INT2, true),

        new SQLColumn('stateId', SQLType.TEXT),
        new SQLColumn('countyId', SQLType.TEXT),

        new SQLColumn('name', SQLType.TEXT),

        new SQLColumn('date', SQLType.BIGINT, true),

        new SQLColumn('confirmed', SQLType.INT, true),
        new SQLColumn('recovered', SQLType.INT, true),
        new SQLColumn('deaths', SQLType.INT, true),
      ],
    });

    axios.get('assets/all_state_rki.json').then((response) => {
      const data = response.data;
      const dataTransformed = [];

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

      this.putAll(dataTransformed);
    });

    axios.get('assets/all_county_rki.json').then((response) => {
      const data = response.data;
      const dataTransformed = [];

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
        this.putAll(subArray);
      }
    });
  }

  async getAllGermany(columns = 'date, confirmed, recovered, deaths') {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.GERMANY} ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getAllState(stateId, columns = 'date, confirmed, recovered, deaths') {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.STATE} AND stateId='${stateId}' ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getAllCounty(countyId, columns = 'date, confirmed, recovered, deaths') {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.COUNTY} AND countyId='${countyId}' ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getAllStates(columns = 'date, stateId, confirmed, recovered, deaths') {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.STATE} ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getAllCountiesOfState(stateId, columns = 'date, countyId, confirmed, recovered, deaths') {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.COUNTY} AND stateId='${stateId}' ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }
}

const RKIStore = new RKISQLStore();
export default RKIStore;
