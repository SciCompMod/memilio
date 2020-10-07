import Datastore from './datastore';
import axios from 'axios';
import _ from 'lodash';

export class Population {
  key;
  name;
  population;
}

/**
 * Interface for simulation data.
 */
class PopulationDatastore extends Datastore {
  constructor() {
    super(
      'population',
      {
        population: '&key, stateKey, name, population',
      },
      (db) => {
        db.population.mapToClass(Population);
        db.population.hook('creating', (pKey, obj, trans) => {
          // set new properties
          obj.name = obj.Bundesland || obj.Landkreis;
          obj.population = obj.EWZ;
          obj.key = obj.stateKey || obj.countyKey;
          if (obj.countyKey) {
            obj.stateKey = parseInt(obj.countyKey / 1000);
          }

          // delete old properties
          delete obj.Bundesland;
          delete obj.Landkreis;
          delete obj.EWZ;
          delete obj.countyKey;
        });
      }
    );
  }

  async populate() {
    const count = await this.db.population.count();
    if (count > 0) {
      return Promise.resolve();
    }

    return new Promise(async (resolve, reject) => {
      let response;
      try {
        response = await axios.get('assets/populations.json');
      } catch (err) {
        // some network error
        reject(err);
        return;
      }

      try {
        await this.db.population.bulkAdd(JSON.parse(JSON.stringify(response.data.states)));
        await this.db.population.bulkAdd(JSON.parse(JSON.stringify(response.data.counties)));

        const germany = _.reduce(response.data.states, (result, value, idx) => result + value.EWZ, 0);

        await this.db.population.add({
          stateKey: -1,
          Bundesland: 'Deutschland',
          EWZ: germany,
        });

        resolve();
      } catch (err) {
        // check error
        if (err.name !== 'BulkError') {
          reject(err);
          return;
        }
      }
    });
  }

  /**
   * Get population for given state or county key
   * @param {number} key state or county key
   * @returns Promise
   */
  async getByKey(key) {
    await this.populate();
    return this.db.population.get(key);
  }

  async getByName(name) {
    await this.populate();
    return this.db.population.where('name').equals(name).first();
  }

  /**
   * Get populations for all counties in state.
   * @param {number} key state key
   * @returns Promise
   */
  async getCountiesByState(key) {
    await this.populate();
    return this.get('population', (e) => e.stateKey === key);
  }

  async clear() {
    this.db.population.clear();
  }
}

const populationstore = new PopulationDatastore();
export default populationstore;
