import Datastore from './datastore';

export class Simulation {
  id;
  created;
  model;
  parameters;
  measures;
  result;
  dataset;
}

/**
 * Interface for simulation data.
 */
class SimulationDatastore extends Datastore {
  constructor() {
    super(
      'simulations',
      {
        simulations: '++, created, model, dataset',
      },
      (db) => {
        db.simulations.mapToClass(Simulation);
        db.simulations.hook('creating', (pKey, obj, trans) => {
          if (!obj.created) {
            obj.created = new Date().getTime();
          }
        });
      }
    );
  }

  /**
   * Get simulation by id
   * @param {number} id Simulation id
   * @returns Promise
   */
  async get(id) {
    return this.db.simulations.get(id);
  }

  /**
   * Get newest simulation result.
   * @returns Promise
   */
  async getLatest() {
    return this.db.simulations.orderBy('created').reverse().first();
  }

  /**
   * Insert a new simulation
   * @param {Simulation} simulation Simulation data
   * @returns Promise
   */
  async add(simulation) {
    return super.add('simulations', simulation);
  }

  /**
   * Update an simulation
   * @param {Simulation} simulation Simulation data
   * @returns Promise
   */
  async put(simulation) {
    return super.put('simulations', simulation);
  }

  /**
   * Delete simulation data
   * @param {number} id Simulation id
   * @returns Promise
   */
  async delete(id) {
    return super.delete('simulations', id);
  }

  /**
   * Delete all simulations stored in database
   * @returns Promise
   */
  async clear() {
    return this.db.simulations.delete();
  }
}

const simulationstore = new SimulationDatastore();
export {simulationstore as SimulationDatastore};
