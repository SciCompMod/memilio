import React, {useState, useEffect} from 'react';
import {Link} from 'react-router-dom';
import axios from 'axios';

export default function Attribution() {
  const [data, setData] = useState({entries: []});

  useEffect(() => {
    document.title = document.title = `Attribution`;
    axios('assets/thirdPartyNotice.json').then((r) => {
      setData({entries: r.data});
    });
  }, []);

  return (
    <div className="attribution">
      <Link tabIndex="1" titel="Zurück zur Hauptseite" to="/">
        Zurück zur Hauptseite
      </Link>
      <h2 className="mt-2">
        <b>Attribution to Third Party Software</b>
      </h2>
      <br />
      <ul>
        {data.entries.map((item) => (
          <li key={item.name + item.version} style={{listStyleType: 'none'}}>
            <h3>
              <b>{item.name + ' ' + item.version}</b>
            </h3>
            <p>Author: {item.author}</p>
            <p>Repository: {item.repository}</p>
            <p>Source: {item.source}</p>
            <p>License: {item.license}</p>
            <p>{item.licenseText}</p>
            <hr className="solid" />
          </li>
        ))}
      </ul>
    </div>
  );
}
