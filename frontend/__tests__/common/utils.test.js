import * as utils from '../../src/common/utils';

it('replace last char', () => {
  expect(utils.replaceLastChar('a,', ';')).toEqual('a;');
});