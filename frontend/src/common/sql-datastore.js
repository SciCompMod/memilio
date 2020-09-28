import DB from './sql-database';
import SQLType from './SQLType';
import {replaceLastChar} from './utils';

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
  /** @type SQLTable */
  table = null;

  /**
   * @param table{SQLTable}
   */
  constructor(table) {
    let createTableQuery = `CREATE TABLE ${table.name} (`;
    for (const column of table.columns) {
      const notNull = column.notNull === true;
      createTableQuery += `${column.name} ${column.type}${notNull ? ' NOT NULL' : ''},`;
    }
    createTableQuery = replaceLastChar(createTableQuery, ');');

    DB.then((db) => db.run(createTableQuery));

    this.table = table;
  }

  async getAll() {
    const query = `SELECT * FROM ${this.table.name};`;
    const result = await DB.then((db) => db.exec(query));
    return this.resultToArray(result);
  }

  /** @param row{any} */
  async put(row) {
    const columns = Object.keys(row);
    const values = Object.values(row);

    let query = `INSERT INTO ${this.table.name} `;

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
    DB.then((db) => db.run(query));
  }

  /**
   * @param {Array<any> | Array<Array<any>>} rows
   * @param {Array<string>} [columns]
   */
  async putAll(rows, columns) {
    let query = `INSERT INTO ${this.table.name} `;

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

    DB.then((db) => db.run(query));
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
