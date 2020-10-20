import {createSlice} from '@reduxjs/toolkit';
import {init} from './app';

const slice = createSlice({
  name: 'models',
  initialState: {},
  reducers: {
    updateModel: (state, action) => {
      const data = {...action.payload};
      if (state[data.model]) {
        let model = state[data.model];

        switch (data.model) {
          case 'secir':
            let parameter = model.parameters.find((p) => p.name === data.parameter);

            if (parameter) {
              parameter[data.group] = parseFloat(data.value);
            }
            break;

          case 'seir':
            model.parameters.find((p) => p.name === data.name).value = data.value;
            break;

          default:
            break;
        }
      }
    },
  },
  extraReducers: {
    [init]: (state, action) => {
      const {models, parameters, ageGroups} = action.payload;

      Object.keys(models).forEach((name) => {
        let model = {
          name,
        };
        if (models[name] instanceof Array) {
          model.parameters = JSON.parse(
            JSON.stringify(
              parameters.filter((p) => models[name].includes(p.name)).map((p) => ({...p, value: p.default}))
            )
          );
        } else if (models[name] instanceof Object) {
          let {ages} = models[name];
          if (ages === 'all') {
            ages = [...ageGroups];
          }
          const params = parameters.filter((p) => models[name].parameters.includes(p.name));

          model.parameters = params.map((p) => {
            return {
              ...p,
              ...ages.reduce((acc, age) => {
                acc[age] = p.default;
                return acc;
              }, {}),
            };
          });
        }
        state[name] = model;
      });
      return state;
    },
  },
});

export const {updateModel} = slice.actions;

export default slice.reducer;
