import axios from 'axios';
import moment from 'moment';
import Datastore from './datastore';
import _ from 'lodash';

const TODAY = moment().startOf('day').valueOf();

class Entry {
  date;
  Confirmed;
  Death;
  Recovered;
}

export class State extends Entry {
  ID_State;
  State;
}

export class County extends Entry {
  ID_County;
  County;
}

export const Tables = {
  STATES: 'states',
  STATES_AGE: 'states_age',
  STATES_GENDER: 'states_gender',
  COUNTIES: 'counties',
  COUNTIES_AGE: 'counties_age',
  COUNTIES_GENDER: 'counties_gender',
};

const Datasets = {
  [Tables.STATES]: 'assets/all_state_rki.json',
  [Tables.STATES_AGE]: 'assets/all_state_age_rki.json',
  [Tables.STATES_GENDER]: 'assets/all_state_gender_rki.json',
  [Tables.COUNTIES]: 'assets/all_county_rki.json',
  [Tables.COUNTIES_AGE]: 'assets/all_county_age_rki.json',
  [Tables.COUNTIES_GENDER]: 'assets/all_county_gender_rki.json',
};

/**
 * Interface for RKI data.
 * All data is stored inside an IndexedDB table.
 */
class RKIDatastore extends Datastore {
  constructor() {
    super(
      'rki',
      {
        latest: '&table, date',
        states: '&[date+ID_State], date, ID_State, State, Confirmed, Death, Recovered',
        states_age: '&[date+ID_State+Age_RKI], date, ID_State, State, Age_RKI, Confirmed, Death, Recovered',
        states_gender: '&[date+ID_State+Gender], date, ID_State, State, Gender, Confirmed, Death, Recovered',
        counties: '&[date+ID_County], date, ID_County, County, ID_State, Confirmed, Death, Recovered',
        counties_age:
          '&[date+ID_County+Age_RKI], date, ID_County, County, Age_RKI, ID_State, Confirmed, Death, Recovered',
        counties_gender:
          '&[date+ID_County+Gender], date, ID_County, County, Gender, ID_State, Confirmed, Death, Recovered',
      },
      (db) => {
        [db.states, db.states_age, db.states_gender].forEach((table) => {
          table.mapToClass(State);
          table.hook('creating', (pkey, obj, trans) => {
            obj.date = obj.Date;
            delete obj.Date;
          });
        });

        [db.counties, db.counties_age, db.counties_gender].forEach((table) => {
          table.mapToClass(County);
          table.hook('creating', (pKey, obj, trans) => {
            obj.ID_State = parseInt(obj.ID_County / 1000);
            obj.date = obj.Date;
            delete obj.Date;
          });
        });
      }
    );
  }

  sanitizeStartEnd(start, end) {
    if (start === null || start === undefined) {
      start = 0;
    } else {
      start = moment(start).startOf('day').valueOf();
    }

    if (end === null || end === undefined) {
      end = Number.MAX_VALUE;
    } else if (end instanceof Date) {
      end = moment(end).endOf('day').valueOf();
    }

    return [start, end];
  }

  populate(tables, progress) {
    if (!Array.isArray(tables)) {
      tables = Array.of(tables);
    }

    if (!(progress instanceof Function)) {
      progress = () => {};
    }

    console.log('populate', tables);

    const {db} = this;

    return new Promise(async (resolve, reject) => {
      const finished = Array(tables.length).fill(false);
      const finish = (idx, value) => {
        finished[idx] = value;

        progress(finished);

        if (finished.every((i) => i !== false)) {
          resolve('finished');
        }
      };

      progress(finished);

      tables
        .map((table) => [table, axios.get(Datasets[table])])
        .forEach(async ([table, request], idx) => {
          let response;
          try {
            response = await request;
          } catch (err) {
            // some network error
            finish(idx, err);
            return;
          }

          try {
            await db[table].bulkAdd(response.data);
          } catch (err) {
            // check error
            if (err.name !== 'BulkError') {
              finish(idx, err);
              return;
            }
          }

          try {
            // try updating timestamp
            await db.latest.put({
              table: table,
              date: TODAY,
            });
          } catch (err) {
            // insert new entry
            await db.latest.add({
              table: table,
              date: TODAY,
            });
          }

          // update germany data
          if (table === Tables.STATES) {
            // delete old germany data
            await this.db.states.where('ID_State').equals(0).delete();

            this.db.transaction('rw', [this.db.states], async (tx) => {
              const dates = await this.db.states.orderBy('date').uniqueKeys();
              dates.forEach(async (date) => {
                const states = await this.db.states.where('date').equals(date).toArray();
                const reduced = _.reduce(
                  states,
                  (result, value, key) => {
                    const {Confirmed, Recovered, Deaths} = value;
                    result.Confirmed += Confirmed;
                    result.Recovered += Recovered;
                    result.Deaths += Deaths;

                    return result;
                  },
                  {
                    Date: date,
                    ID_State: 0,
                    State: 'Deutschland',
                    Confirmed: 0,
                    Recovered: 0,
                    Deaths: 0,
                  }
                );

                this.db.states.add(reduced);
              });
            });
          }

          finish(idx, true);
        });
    });
  }

  async get(table, options, filter) {
    // check if table is up to date
    let latest = undefined;
    try {
      latest = await this.db.latest.get(table);
    } catch (err) {}

    if (latest === undefined || moment(latest.date).isBefore(TODAY)) {
      // no current data for table, so download latest data
      await this.populate(table);
    }

    const {start, end, sort, where} = options;
    const [s, e] = this.sanitizeStartEnd(start, end);

    let query = this.db[table];
    if (where !== undefined && Array.isArray(where) && where.length === 2) {
      query = where[1](query.where(where[0])).and((entry) => entry.date >= s && entry.date <= e);
    } else {
      query = query.where('date').between(s, e, true, true);
    }

    if (filter !== undefined) {
      query = query.and(filter);
    }

    if (sort) {
      return query.sortBy(sort);
    }

    return query.toArray();
  }

  /**
   * Retrives case data for given state for specified interval.
   *
   * @param {number} stateID ID of state to be retrived
   * @param {number | string | Date} start start date of interval to be retrived (inclusive)
   * @param {number | string | Date} end end date of interval to be retrived (inclusive)
   *
   * @returns Promise
   */
  getState(stateID, options = {}) {
    const {age, gender} = options;

    if (stateID === undefined) {
      Promise.reject('No state id provided!');
    }

    if (age !== undefined && gender !== undefined) {
      Promise.reject('Only age or gender can be selected!');
    }

    const where = ['ID_State', (wc) => wc.equals(stateID)];

    let promise;

    if (age !== undefined) {
      promise = this.get(Tables.STATES_AGE, {...options, where}, (e) => {
        return age instanceof String ? e.Age_RKI === age : true;
      });
    } else if (gender !== undefined) {
      promise = this.get(Tables.STATES_GENDER, {...options, where}, (e) => {
        return gender instanceof String ? e.Gender === gender : true;
      });
    } else {
      promise = this.get(Tables.STATES, {...options, where});
    }

    return promise;
  }

  /**
   * Retrives case data for given county for specified interval.
   *
   * @param {number} countyID ID of county to be retrived
   * @param {number | string | Date} start start date of interval to be retrived (inclusive)
   * @param {number | string | Date} end end date of interval to be retrived (inclusive)
   *
   * @returns Promise
   */
  getCounty(countyID, options = {}) {
    const {age, gender} = options;

    if (countyID === undefined) {
      Promise.reject('No county id provided!');
    }

    if (age !== undefined && gender !== undefined) {
      Promise.reject('Only age or gender can be selected!');
    }

    const where = ['ID_County', (wc) => wc.equals(countyID)];

    let promise;

    if (age !== undefined) {
      promise = this.get(Tables.COUNTIES_AGE, {...options, where}, (e) => {
        return age instanceof String ? e.Age_RKI === age : true;
      });
    } else if (gender !== undefined) {
      promise = this.get(Tables.COUNTIES_GENDER, {...options, where}, (e) => {
        return gender instanceof String ? e.Gender === gender : true;
      });
    } else {
      promise = this.get(Tables.COUNTIES, {...options, where});
    }

    return promise;
  }

  /**
   * Retrives case data for all counties in given state for specified interval.
   *
   * @param {number} stateID ID of state to be retrived
   * @param {number | string | Date} start start date of interval to be retrived (inclusive)
   * @param {number | string | Date} end end date of interval to be retrived (inclusive)
   *
   * @returns Promise
   */
  getCountiesByState(stateID, options = {}) {
    const {age, gender} = options;

    if (stateID === undefined || stateID === null) {
      Promise.reject('No state id provided!');
    }

    if (age !== undefined && gender !== undefined) {
      Promise.reject('Only age or gender can be selected!');
    }

    let promise;

    const where = ['ID_State', (wc) => wc.equals(stateID)];

    if (age !== undefined) {
      promise = this.get(Tables.COUNTIES_AGE, {...options, where}, (e) => {
        return age instanceof String ? e.Age_RKI === age : true;
      });
    } else if (gender !== undefined) {
      promise = this.get(Tables.COUNTIES_GENDER, {...options, where}, (e) => {
        return gender instanceof String ? e.Gender === gender : true;
      });
    } else {
      promise = this.get(Tables.COUNTIES, {...options, where});
    }

    return promise;
  }

  async getTimeBounds(table) {
    const start = await this.db[table].orderBy('date').first();
    const end = await this.db[table].orderBy('date').last();

    if (start.date < new Date('2020-02-24').getTime()) {
      start.date = moment('2020-02-24').startOf('day').valueOf();
    }

    return {start: start.date, end: end.date};
  }

  /**
   * Clear all data in tables, but does not delete the tables.
   *
   * @returns Promise
   */
  async clear() {
    return Promise.all(
      [
        this.db.latest,
        this.db.states,
        this.db.states_age,
        this.db.states_gender,
        this.db.counties,
        this.db.counties_age,
        this.db.counties_gender,
      ].map((table) => table.toCollection().delete())
    );
  }

  /**
   * Delete database. It cannot be used afterwards anymore!
   */
  async delete() {
    return this.db.delete();
  }
}

const rkidatastore = new RKIDatastore();
export {rkidatastore as RKIDatastore};
