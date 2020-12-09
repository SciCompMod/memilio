import SQLType from '../../../../common/datastore/sql/SQLType';

it('DECIMAL', () => {
  expect(SQLType.DECIMAL(3, 5)).toEqual('DECIMAL(3,5)');
  expect(SQLType.DECIMAL(1, 2)).toEqual('DECIMAL(1,2)');
  expect(SQLType.DECIMAL(1, 0)).toEqual('DECIMAL(1,0)');
});

it('CHARACTER', () => {
  expect(SQLType.CHARACTER(1)).toEqual('CHARACTER(1)');
  expect(SQLType.CHARACTER(2)).toEqual('CHARACTER(2)');
  expect(SQLType.CHARACTER(128)).toEqual('CHARACTER(128)');
});

it('VARCHAR', () => {
  expect(SQLType.VARCHAR(1)).toEqual('VARCHAR(1)');
  expect(SQLType.VARCHAR(2)).toEqual('VARCHAR(2)');
  expect(SQLType.VARCHAR(128)).toEqual('VARCHAR(128)');
});

it('VARYING_CHARACTER', () => {
  expect(SQLType.VARYING_CHARACTER(1)).toEqual('VARYING CHARACTER(1)');
  expect(SQLType.VARYING_CHARACTER(2)).toEqual('VARYING CHARACTER(2)');
  expect(SQLType.VARYING_CHARACTER(128)).toEqual('VARYING CHARACTER(128)');
});

it('VARYING_CHARACTER', () => {
  expect(SQLType.NCHAR(1)).toEqual('NCHAR(1)');
  expect(SQLType.NCHAR(2)).toEqual('NCHAR(2)');
  expect(SQLType.NCHAR(128)).toEqual('NCHAR(128)');
});

it('VARYING_CHARACTER', () => {
  expect(SQLType.NATIVE_CHARACTER(1)).toEqual('NATIVE CHARACTER(1)');
  expect(SQLType.NATIVE_CHARACTER(2)).toEqual('NATIVE CHARACTER(2)');
  expect(SQLType.NATIVE_CHARACTER(128)).toEqual('NATIVE CHARACTER(128)');
});

it('VARYING_CHARACTER', () => {
  expect(SQLType.NVARCHAR(1)).toEqual('NVARCHAR(1)');
  expect(SQLType.NVARCHAR(2)).toEqual('NVARCHAR(2)');
  expect(SQLType.NVARCHAR(128)).toEqual('NVARCHAR(128)');
});
