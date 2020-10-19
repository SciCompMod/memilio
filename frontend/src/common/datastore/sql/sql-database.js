import initSqlJs from 'sql.js';

const promise = new Promise((resolve, reject) => {
  initSqlJs({
    locateFile: (file) => `wasm/${file}`,
  })
    .then((SQL) => {
      const db = new SQL.Database();
      resolve(db);
    })
    .catch(reject);
});

export default promise;
