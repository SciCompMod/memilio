import Dexie from 'dexie';

/**
 * Interface for interacting with IndexedDB database.
 */
export default class Datastore {
  /**
   *
   * @param {String} name Name used for IndexedDB database
   * @param {Object} tables Table schema definitions (see https://dexie.org/docs/Version/Version.stores() for informations about schema definition)
   * @param {Function} init Function for additional initialization logic
   * @param {Function} upgrade Function for database upgrade
   * @param {Number} version Newest database version number
   */
  constructor(name, tables, init = () => {}, upgrade = []) {
    this.db = new Dexie(name);
    this.db.delete();
    if (!Array.isArray(tables)) {
      tables = [tables];
    }

    tables.forEach((t, idx) => {
      let version = this.db.version(idx + 1).stores(t);
      if (idx > 0 && upgrade.length > idx && upgrade[idx] !== null) {
        version.upgrade(upgrade[idx]);
      }
    });

    init(this.db);

    this.db.open();
  }

  /**
   * Basic method to retrieve values from IndexedDB table
   *
   * @param {String} table IndexedDB table to retrieve data from
   * @param {Function} filter Filter Function
   * @param {String} sort Column name to sort by. For descending sort prepend name with '-'.
   */
  async get(table, filter, sort) {
    let collection = this.db[table].filter(filter || (() => true));

    if (sort) {
      if (sort.startsWith('-')) {
        collection = collection.reverse();
        sort = sort.substring(1);
      }
      collection = collection.sortBy(sort);
    } else {
      collection = collection.toArray();
    }

    return collection;
  }

  /**
   * Basic method to add values to IndexedDB table
   *
   * @param {String} table IndexedDB table to retrieve data from
   * @param {Oject} data Data to store
   */
  async add(table, entry) {
    return this.db[table].add(entry);
  }

  /**
   * Basic method update values inside IndexedDB table
   *
   * @param {String} table IndexedDB table to retrieve data from
   * @param {Object} data Data entry to update
   */
  async put(table, data) {
    return this.db[table].put(data);
  }

  /**
   * Basic method to delete values from IndexedDB table
   *
   * @param {String} table IndexedDB table to retrieve data from
   * @param {Function} primaryKey Primary Key value of entry to delete
   */
  async delete(table, primaryKey) {
    return this.db[table].delete(primaryKey);
  }

  async clear() {}
}
