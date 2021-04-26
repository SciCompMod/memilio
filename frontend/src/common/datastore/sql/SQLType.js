/**
 * @callback SQLTypeCallback
 * @return string
 */

/**
 * @readonly
 * @enum {string | SQLTypeCallback}
 */
const SQLType = {
  NULL: 'NULL',

  INT: 'INT',
  INTEGER: 'INTEGER',
  TINYINT: 'TINYINT',
  SMALLINT: 'SMALLINT',
  MEDIUMINT: 'MEDIUMINT',
  BIGINT: 'BIGINT',
  UNSIGNED_BIG_INT: 'UNSIGNED BIG INT',
  INT2: 'INT2',
  INT8: 'INT8',

  REAL: 'REAL',
  DOUBLE: 'DOUBLE',
  DOUBLE_PRECISION: 'DOUBLE PRECISION',
  FLOAT: 'FLOAT',

  NUMERIC: 'NUMERIC',

  DECIMAL: (a, b) => `DECIMAL(${a},${b})`,
  BOOLEAN: 'BOOLEAN',
  DATE: 'DATE',
  DATETIME: 'DATETIME',

  CHARACTER: (size) => `CHARACTER(${size})`,
  VARCHAR: (size) => `VARCHAR(${size})`,
  VARYING_CHARACTER: (size) => `VARYING CHARACTER(${size})`,
  NCHAR: (size) => `NCHAR(${size})`,
  NATIVE_CHARACTER: (size) => `NATIVE CHARACTER(${size})`,
  NVARCHAR: (size) => `NVARCHAR(${size})`,
  TEXT: 'TEXT',
  CLOB: 'CLOB',

  BLOB: 'BLOB',
};

export default SQLType;
