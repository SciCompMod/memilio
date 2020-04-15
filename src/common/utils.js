
const timeParser = d3.timeParse("%d.%m.%Y");

function parseTime(date_string) {
    return timeParser(date_string);
}

export {
    parseTime
}