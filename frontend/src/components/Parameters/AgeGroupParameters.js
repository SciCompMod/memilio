import React, {Component} from 'react';
import {withTranslation} from 'react-i18next';

import {kaReducer, Table} from 'ka-table';
import {DataType, EditingMode} from 'ka-table/enums';

import 'ka-table/style.scss';
import './AgeGroupParameters.scss';

const tablePropsInit = {
  columns: [],
  data: [],
  editableCells: [],
  editingMode: EditingMode.Cell,
  rowKeyField: 'name',
};

class AgeGroupParameters extends Component {
  constructor(props) {
    super(props);

    const {parameters, ageGroups, t} = this.props;
    const tableProps = JSON.parse(JSON.stringify(tablePropsInit));

    tableProps.columns = [
      {
        dataType: DataType.String,
        key: 'label',
        title: t('parameters.title'),
        format: (value) => t(`${value}.label`),
        isEditable: false,
      },
    ];

    ageGroups.forEach((group) => {
      tableProps.columns.push({
        dataType: DataType.String,
        key: group,
        title: group,
        validation: (value, rowData) => {
          if (value < rowData.min || rowData.max < value) {
            // TODO use translation service
            return `Value must be between ${rowData.min} and ${rowData.max}`;
          }
        },
      });
    });

    tableProps.data = [...parameters];

    this.state = {
      tableProps,
    };
  }

  dispatch(action) {
    if (action.type === 'UpdateCellValue') {
      this.props.onChange({
        parameter: action.rowKeyValue,
        group: action.columnKey,
        value: action.value,
      });
    }

    this.setState((state) => {
      if (state.tableProps.editableCells.length > 0) {
        state.tableProps.editableCells = [];
      }

      return {
        tableProps: kaReducer(state.tableProps, action),
      };
    });
  }

  render() {
    return (
      <div className="age-group-parameters">
        <Table
          {...this.state.tableProps}
          childAttributes={{
            table: {
              className: 'age-group-table',
            },
            cell: {
              className: 'age-group-cell',
            },
          }}
          dispatch={this.dispatch.bind(this)}
        />
      </div>
    );
  }
}

export default withTranslation()(AgeGroupParameters);
