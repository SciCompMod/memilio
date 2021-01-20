import React, {useEffect} from 'react';
import {Link} from 'react-router-dom';

export default function About() {
  useEffect(() => {
    document.title = document.title = `Über die Webseite`;
  });

  return (
    <div className="about">
      <Link tabIndex="1" titel="Zurück zur Hauptseite" to="/">
        Zurück zur Hauptseite
      </Link>
      <h1 className="mt-2">Über die Webseite</h1>
      <p>
        Die Corona-Pandemie stellt nicht nur jeden einzelnen Menschen vor viele Fragen. Auch die Politik muss
        Entscheidungen über sinnvolle und vertretbare Gegenmaßnahmen und Eindämmungsversuche treffen. Für politische und
        wirtschaftliche Entscheidungen werden fundierte Informationen aus Bereichen der Virologie und Epidemiologie
        benötigt. Darüber hinaus kann aber auch die Expertise in den Bereichen der Mathematik und Informatik, speziell
        der numerischen Simulation einen wichtigen Beitrag leisten.
      </p>
      <p>
        Das Institut für Softwaretechnologie entwickelt in Zusammenarbeit mit dem Helmholtz-Zentrum für
        Infektionsforschung ein umfassendes Softwarepaket, welches das COVID19-Infektionsgeschehen per Simulation
        darstellt. Neben einer hochaufgelösten Modellierung der Infektionsausbreitung steht im Mittelpunkt die Wirkung
        der Gegenmaßnahmen zu bewerten. Teil des Softwarepaktes ist ein Online-Tool, welches der Öffentlichkeit zu
        Verfügung gestellt wird. Das Tool soll zum einen der Politik dienen, um sich bei Entscheidungen zu weiteren
        Maßnahmen auf fundierten wissenschaftlichen Daten und Hochrechnungen zu stützen. Zum anderen soll das intuitive
        User-Interface helfen, den Verlauf der Infektionsketten auch der Bevölkerung erkennbar zu machen.
      </p>
      <p>
        In der hier veröffentlichen Visualisierung werden die Inzidenzwerte der dem Robert Koch-Institut gemeldeten
        Fälle pro 7 Tage und 100.000 Einwohner und die aktuelle Reproduktionszahl in den einzelnen Landkreisen
        angezeigt. Die Reproduktionszahl gibt an, wie viele Menschen unter den aktuellen Maßnahmen von einer infektiösen
        Person durchschnittlich angesteckt werden. Eine Untersuchung dieser Größe auf Landkreisebene bekommt durch das
        Auftreten einer Virusvariante mit deutlich erhöhter Übertragbarkeit eine neue Bedeutung, denn eine auffällig
        erhöhte lokale Reproduktionszahl könnte auf einen Einfluss dieser Virusvariante hindeuten. Deren
        Reproduktionszahl ist nach Untersuchungen des Imperial College London (
        <a
          href="https://www.medrxiv.org/content/10.1101/2020.12.30.20249034v2"
          target="_blank"
          rel="noopener noreferrer"
        >
          Link
        </a>
        ) deutlich höher und sie kann daher zu einer beschleunigten lokalen Ausbreitung führen. Die Visualisierung
        basiert auf den dem Robert Koch-Institut gemeldeten Fällen. Die <b>absolute</b> Visualisierung gibt die
        berechnete aktuelle Reproduktionszahl pro Landkreis an. Die <b>relative</b> Visualisierung gibt den Quotienten
        der Reproduktionszahl des Landkreises durch die Reproduktionszahl für ganz Deutschland aus. Zu beachten ist,
        dass die Berechnung der Reproduktionszahl auf Landkreisebene größeren statistischen Schwankungen unterliegen
        kann. Dies ist insbesondere der Fall, wenn die die Fallzahlen, die der Berechnung zu Grunde liegen, sehr niedrig
        sind.
      </p>
      <p>
        <b>Disclaimer:</b> This website and its contents herein, including all data and analysis are provided to the
        public strictly for educational and academic research purposes. Reliance on the Website for medical guidance or
        use of the Website in commerce is strictly prohibited.
      </p>
      <h2>Ansprechpartner Helmholtz-Zentrum für Infektionsforschung</h2>
      <p>
        Prof. Dr. Michael Meyer-Hermann
        <br />
        Dr. Sebastian Binder
        <br />
        <a
          href="https://www.helmholtz-hzi.de/de/forschung/forschungsschwerpunkte/immunantwort-und-interventionen/system-immunologie/m-meyer-hermann/"
          target="_blank"
          rel="noopener noreferrer"
        >
          Kontakt
        </a>
      </p>
      <h2>Ansprechpartner Deutsches Zentrum für Luft- und Raumfahrt</h2>
      Dr. Martin J. Kühn
      <br />
      Dr. Margrit Klitz
      <br />
      <a href="mailto:contact-hpc-against-corona@dlr.de">Kontakt</a>
    </div>
  );
}
