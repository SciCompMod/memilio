import * as utils from '../../src/common/utils';
import {lastElement, roundToUTCNoon} from '../../src/common/utils';

it('roundToUTCMidnight', () => {
  const input0 = Date.UTC(2000, 6, 15, 0, 0, 0, 0);
  const expected0 = Date.UTC(2000, 6, 15, 12, 0, 0, 0);
  expect(roundToUTCNoon(input0)).toEqual(expected0);

  const input1 = Date.UTC(2000, 6, 15, 12, 0, 0, 0);
  const expected1 = Date.UTC(2000, 6, 15, 12, 0, 0, 0);
  expect(roundToUTCNoon(input1)).toEqual(expected1);

  const input2 = Date.UTC(2000, 6, 15, 23, 59, 59, 999);
  const expected2 = Date.UTC(2000, 6, 15, 12, 0, 0, 0);
  expect(roundToUTCNoon(input2)).toEqual(expected2);
});

it('filterJSObject', () => {
  // TODO
});

it('stateIdFromCountyId', () => {
  expect(utils.stateIdFromCountyId('01234')).toEqual('01');
  expect(utils.stateIdFromCountyId('16234')).toEqual('16');

  expect(() => utils.stateIdFromCountyId('')).toThrowError();
  expect(() => utils.stateIdFromCountyId('9')).toThrowError();
  expect(() => utils.stateIdFromCountyId('16')).toThrowError();
  expect(() => utils.stateIdFromCountyId('6234')).toThrowError();
  expect(() => utils.stateIdFromCountyId('a1234')).toThrowError();
  expect(() => utils.stateIdFromCountyId('1234a')).toThrowError();
});

it('isStateId', () => {
  expect(utils.isStateId('00')).toBeTruthy();
  expect(utils.isStateId('01')).toBeTruthy();
  expect(utils.isStateId('02')).toBeTruthy();
  expect(utils.isStateId('12')).toBeTruthy();

  expect(utils.isStateId('')).toBeFalsy();
  expect(utils.isStateId('0')).toBeFalsy();
  expect(utils.isStateId('-1')).toBeFalsy();
  expect(utils.isStateId('+1')).toBeFalsy();
  expect(utils.isStateId('1+')).toBeFalsy();
  expect(utils.isStateId('-10')).toBeFalsy();
  expect(utils.isStateId('+10')).toBeFalsy();
  expect(utils.isStateId('000')).toBeFalsy();
  expect(utils.isCountyId('aa')).toBeFalsy();
  expect(utils.isCountyId('x0')).toBeFalsy();
  expect(utils.isCountyId('d0')).toBeFalsy();
  expect(utils.isCountyId('0d')).toBeFalsy();
  expect(utils.isCountyId('0b')).toBeFalsy();
  expect(utils.isCountyId('0x00')).toBeFalsy();
  expect(utils.isCountyId('0d00')).toBeFalsy();
  expect(utils.isCountyId('0b00')).toBeFalsy();
});

it('isCountyId', () => {
  expect(utils.isCountyId('00000')).toBeTruthy();
  expect(utils.isCountyId('12345')).toBeTruthy();
  expect(utils.isCountyId('99999')).toBeTruthy();

  expect(utils.isCountyId('')).toBeFalsy();
  expect(utils.isCountyId('0')).toBeFalsy();
  expect(utils.isCountyId('00')).toBeFalsy();
  expect(utils.isCountyId('000')).toBeFalsy();
  expect(utils.isCountyId('0000')).toBeFalsy();
  expect(utils.isCountyId('aaaaa')).toBeFalsy();
  expect(utils.isCountyId('0000a')).toBeFalsy();
  expect(utils.isCountyId('-1000')).toBeFalsy();
  expect(utils.isCountyId('+1000')).toBeFalsy();
  expect(utils.isCountyId('-10000')).toBeFalsy();
  expect(utils.isCountyId('+10000')).toBeFalsy();
  expect(utils.isCountyId('00a00')).toBeFalsy();
  expect(utils.isCountyId('0x000')).toBeFalsy();
  expect(utils.isCountyId('0d000')).toBeFalsy();
  expect(utils.isCountyId('0b000')).toBeFalsy();
  expect(utils.isCountyId('0x00000')).toBeFalsy();
  expect(utils.isCountyId('0d00000')).toBeFalsy();
  expect(utils.isCountyId('0b00000')).toBeFalsy();
  expect(utils.isCountyId('000000')).toBeFalsy();
});

it('lastElement', () => {
  expect(lastElement([0])).toEqual(0);
  expect(lastElement([0, 0])).toEqual(0);
  expect(lastElement([0, 1])).toEqual(1);
  expect(() => lastElement([])).toThrow(RangeError);
});

it('replaceLastChar', () => {
  expect(utils.replaceLastChar('a,', ';')).toEqual('a;');
  expect(utils.replaceLastChar(',', ';')).toEqual(';');
  expect(() => utils.replaceLastChar('', ';')).toThrow(RangeError);
});
