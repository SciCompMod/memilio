import SQLDatastore, {SQLColumn} from './sql-datastore';
import SQLType from './SQLType';
import {stateIdFromCountyId} from '../../utils';

class SECIRSQLStore extends SQLDatastore {
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
      name: 'secir',
      columns: [
        new SQLColumn('level', SQLType.INT2, true),

        new SQLColumn('stateId', SQLType.INT),
        new SQLColumn('countyId', SQLType.INT),

        new SQLColumn('name', SQLType.TEXT),

        new SQLColumn('date', SQLType.BIGINT, true),

        new SQLColumn('susceptible', SQLType.INT, true),
        new SQLColumn('exposed', SQLType.INT, true),
        new SQLColumn('infected', SQLType.INT, true),
        new SQLColumn('recovered', SQLType.INT, true),
        new SQLColumn('car', SQLType.INT, true),
        new SQLColumn('hospitalized', SQLType.INT, true),
        new SQLColumn('icu', SQLType.INT, true),
        new SQLColumn('dead', SQLType.INT, true),
      ],
    };

    this.create(this.table);
  }

  async populate(level, id, name, data) {
    await this.clear('secir');

    let l;
    let stateId = null;
    let countyId = null;
    switch (level) {
      case 'germany':
        l = this.Level.GERMANY;
        break;
      case 'states':
        l = this.Level.STATE;
        stateId = id;
        break;
      case 'counties':
        l = this.Level.COUNTY;
        stateId = stateIdFromCountyId(id);
        countyId = id;
        break;
      default:
        break;
    }

    let dataTransformed = [];
    for (const datum of data) {
      dataTransformed.push([
        l,
        stateId,
        countyId,
        name,
        datum.date,
        datum.S,
        datum.E,
        datum.I,
        datum.R,
        datum.nb_car,
        datum.nb_hosp,
        datum.nb_icu,
        datum.nb_dead,
      ]);
    }

    await this.putAll('secir', dataTransformed);
  }

  async getAllGermany() {
    // const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.GERMANY} ORDER BY date`;
    const query = `SELECT 
      date, 
      SUM(susceptible) as susceptible, 
      SUM(exposed) as exposed, 
      SUM(infected) as infected 
      SUM(recovered) as recovered 
      SUM(car) as car
      SUM(hospitalized) as hospitalized 
      SUM(icu) as icu 
      SUM(dead) as dead 
    FROM ${this.table.name}
    WHERE level=${this.Level.STATE}
    GROUP BY date
    ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getAllState(
    stateId,
    columns = 'date, susceptible, exposed, infected, recovered, car, hospitalized, icu, dead'
  ) {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.STATE} AND stateId='${stateId}' ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getAllStateInRange(
    stateId,
    range,
    columns = 'date, susceptible, exposed, infected, recovered, car, hospitalized, icu, dead'
  ) {
    const query = `SELECT 
      ${columns} 
    FROM ${this.table.name} 
    WHERE level=${this.Level.STATE} 
      AND stateId='${stateId}' 
      AND (date BETWEEN ${range.start} AND ${range.end}) 
    ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getAllCounty(
    countyId,
    columns = 'date, susceptible, exposed, infected, recovered, car, hospitalized, icu, dead'
  ) {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.COUNTY} AND countyId='${countyId}' ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getAllStates(
    columns = 'date, stateId, susceptible, exposed, infected, recovered, car, hospitalized, icu, dead'
  ) {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.STATE} ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getAllStatesInRange(
    range,
    columns = 'date, stateId, susceptible, exposed, infected, recovered, car, hospitalized, icu, dead'
  ) {
    const query = `SELECT 
      ${columns} 
    FROM ${this.table.name} 
    WHERE level=${this.Level.STATE}  
      AND (date BETWEEN ${range.start} AND ${range.end}) 
    ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getAllCountiesOfState(
    stateId,
    columns = 'date, countyId, susceptible, exposed, infected, recovered, car, hospitalized, icu, dead'
  ) {
    const query = `SELECT ${columns} FROM ${this.table.name} WHERE level=${this.Level.COUNTY} AND stateId='${stateId}' ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getAllCountiesInRange(
    range,
    columns = 'date, countyId, susceptible, exposed, infected, recovered, car, hospitalized, icu, dead'
  ) {
    const query = `SELECT 
      ${columns} 
    FROM ${this.table.name} 
    WHERE level=${this.Level.COUNTY}  
      AND (date BETWEEN ${range.start} AND ${range.end}) 
    ORDER BY date`;
    return this.resultToArray(await this.execQuery(query));
  }

  async getTimeBounds() {
    const query = `SELECT MIN(date) as start, MAX(date) as end FROM ${this.table.name} WHERE level=${this.Level.STATE}`;
    return this.resultToArray(await this.execQuery(query))[0];
  }
}

const SECIRStore = new SECIRSQLStore();
export default SECIRStore;
