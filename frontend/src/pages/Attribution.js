import React, {useState, useEffect} from 'react';
import {Link, useRouteMatch} from 'react-router-dom';
import axios from 'axios';

import {getParentRoute} from '~/common/utils';

export default function Attribution() {
  const {path} = useRouteMatch();
  const [data, setData] = useState({entries: []});

  useEffect(() => {
    document.title = `Attribution`;
    axios('assets/thirdPartyNotice.json').then((r) => {
      setData({entries: r.data});
    });
  }, []);

  return (
    <div className="attribution">
      <Link tabIndex="1" titel="Zurück zur Hauptseite" to={getParentRoute(path)}>
        Zurück zur Hauptseite
      </Link>
      <h2 className="mt-2">
        <b>Attribution to Third Party Software</b>
      </h2>
      <br />
      <p>
        Several fantastic pieces of free and open-source software have really helped us to get the website to where it
        is today! We are grateful, that this is possible in today's world and collaboration around the globe is stronger
        than ever before.
      </p>
      <hr className="solid" />
      <ul>
        {data.entries.map((item) => (
          <li key={item.name + item.version} style={{listStyleType: 'none'}}>
            <h3>
              <b>
                <u>{item.name + ' ' + item.version}</u>
              </b>
            </h3>
            <p>
              <b>Author:</b> {item.author}
            </p>
            <p>
              <b>Repository:</b> {item.repository}
            </p>
            <p>
              <b>Source:</b> {item.source}
            </p>
            <p>
              <b>License:</b> {item.license}
            </p>
            {item.licenseText?.split('\n\n').map((paragraph) => (
              <p>{paragraph}</p>
            ))}
            <hr className="solid" />
          </li>
        ))}
      </ul>
    </div>
  );
}
