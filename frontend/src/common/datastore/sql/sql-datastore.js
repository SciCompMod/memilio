import DB from './sql-database';
import {replaceLastChar} from '../../utils';

export class SQLColumn {
  /** @type string */
  name;

  /** @type SQLType */
  type;

  /** @type boolean */
  notNull = false;

  /**
   * @param name{string}
   * @param type{SQLType}
   * @param [notNull]{boolean}
   */
  constructor(name, type, notNull = false) {
    this.name = name;
    this.type = type;
    this.notNull = notNull;
  }
}

export class SQLTable {
  /** @type string | null */
  name = null;

  /** @type Array<SQLColumn> */
  columns = [];
}

/**
 * An abstraction around an in memory SQL database.
 */
export default class SQLDatastore {
  /**
   * Abstract method to load data into datastore
   */
  async populate() {
    return false;
  }

  /**
   * @param table{SQLTable}
   */
  async create(table) {
    let createTableQuery = `CREATE TABLE ${table.name} (`;
    for (const column of table.columns) {
      const notNull = column.notNull === true;
      createTableQuery += `${column.name} ${column.type}${notNull ? ' NOT NULL' : ''},`;
    }
    createTableQuery = replaceLastChar(createTableQuery, ');');

    await this.runQuery(createTableQuery);
  }

  /**
   * Returns all entries of the given table.
   * @param table{string} The name of the table.
   * @return Array<any> Each array entry is an object representing a row,
   *                    with keys corresponding to column names and values to the values.
   */
  async getAll(table) {
    const query = `SELECT * FROM ${table};`;
    const result = await DB.then((db) => db.exec(query));
    return this.resultToArray(result);
  }

  /**
   * @param table{string} The name of the table to put the row in
   * @param row{any} A JavaScript object where the keys are the column names and the values are the values.
   */
  async put(table, row) {
    const columns = Object.keys(row);
    const values = Object.values(row);

    let query = `INSERT INTO ${table} `;

    query += '(';
    for (const column of columns) {
      query += `${column},`;
    }
    query = replaceLastChar(query, ') ');

    query += 'VALUES (';
    for (const value of values) {
      query += `'${value}',`;
    }
    query = replaceLastChar(query, ');');
    await this.runQuery(query);
  }

  /**
   * Inserts multiple rows into the given table. The columns parameter has to be specified only if you don't insert all
   * columns in the order they have in the table.
   *
   * @param table{string} The name of the table.
   * @param {Array<any> | Array<Array<any>>} rows
   * @param {Array<string>} [columns]
   */
  async putAll(table, rows, columns) {
    let query = `INSERT INTO ${table} `;

    if (columns) {
      query += '(';
      for (const column of columns) {
        query += `${column},`;
      }
      query = replaceLastChar(query, ') ');
    }

    query += 'VALUES ';
    for (const row of rows) {
      query += '(';
      for (const value of Object.values(row)) {
        query += `'${value}',`;
      }
      query = replaceLastChar(query, '),');
    }
    query = replaceLastChar(query, ';');

    await this.runQuery(query);
  }

  /**
   * Deletes all entries of the table and frees used memory.
   * @param table{string} The name of the table.
   * @return {Promise<void>}
   */
  async clear(table) {
    await this.runQuery(`DELETE FROM ${table}; VACUUM`);
  }

  /**
   * Runs the given query and ignores the result.
   * @param query{string}
   * @return {Promise<void>}
   */
  async runQuery(query) {
    DB.then((db) => db.run(query));
  }

  /**
   * Runs the given query and returns the result.
   * @param query{string}
   * @return {Promise<Array<{columns: Array<string>, values: Array<Array<any>>}>>}
   */
  async execQuery(query) {
    return await DB.then((db) => db.exec(query));
  }

  /**
   * @protected
   *
   * Converts the result of a query to an array of rows, where rows are objects with their key representing column names
   * and their values representing entries.
   *
   * @param {Array<{columns: Array<string>, values: Array<Array<any>>}>} result
   * @return Array<any>
   */
  resultToArray(result) {
    if (result.length === 0) return [];

    const columns = result[0].columns;
    const rows = result[0].values;
    const array = [];

    for (const row of rows) {
      const obj = {};
      for (let i = 0; i < columns.length; ++i) {
        obj[columns[i]] = row[i];
      }
      array.push(obj);
    }

    return array;
  }
}
