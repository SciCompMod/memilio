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

  async getAll(table) {
    const query = `SELECT * FROM ${table};`;
    const result = await DB.then((db) => db.exec(query));
    return this.resultToArray(result);
  }

  /** @param row{any} */
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

  async clear(table) {
    await this.runQuery(`DELETE FROM ${table}`);
  }

  async runQuery(query) {
    DB.then((db) => db.run(query));
  }

  async execQuery(query) {
    return await DB.then((db) => db.exec(query));
  }

  /**
   * @protected
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
