import React from 'react';
import {Link} from 'react-router-dom';

export default function About() {
  return (
    <div className="about">
      <Link to="/">Zurück</Link>
      <h1 className="mt-2">Über die Webseite</h1>
      <p>
        Die Corona-Pandemie stellt nicht nur jeden einzelnen Menschen vor viele Fragen. Auch die Politik muss
        Entscheidungen über sinnvolle und vertretbare Gegenmaßnahmen und Eindämmungsversuche treffen. Für politische und
        wirtschaftliche Entscheidungen werden fundierte Informationen aus Bereichen der Virologie und Epidemiologie
        benötigt. Darüber hinaus kann aber auch die Expertise in den Bereichen der Mathematik und Informatik, speziell
        der numerischen Simulation einen wichtigen Beitrag leisten.
      </p>
      <p>
        Das Institut für Softwaretechnologie entwickelt in Zusammenarbeit mit dem Helmholtzzentrum für
        Infektionsforschung ein umfassendes Softwarepaket, welches das COVID19-Infektionsgeschehen per Simulation
        darstellt. Neben einer hochaufgelösten Modellierung der Infektionsausbreitung steht im Mittelpunkt die Wirkung
        der Gegenmaßnahmen zu bewerten. Teil des Softwarepaktes ist ein Online-Tool, welches der Öffentlichkeit zu
        Verfügung gestellt wird. Das Tool soll zum einen der Politik dienen, um sich bei Entscheidungen zu weiteren
        Maßnahmen auf fundierten wissenschaftlichen Daten und Hochrechnungen zu stützen. Zum anderen soll das intuitive
        User-Interface helfen, den Verlauf der Infektionsketten auch der Bevölkerung erkennbar zu machen.
      </p>
      <p>
        In der hier veröffentlichen Visualisierung wird die aktuelle Reproduktionszahl in den einzelnen Landkreisen
        angezeigt. Die Reproduktionszahl gibt an, wie viele Menschen unter den aktuellen Maßnamen von einer infektiösen
        Person durchschnittlich angesteckt werden. Eine Untersuchung dieser Größe auf Landkreisebene bekommt durch das
        Auftreten einer Virusvariante mit deutlich erhöhter Übertragbarkeit eine neue Bedeutung, denn eine auffällig
        erhöhte lokale Reproduktionszahl könnte auf einen Einfluss dieser Virusvariante hindeuten. Deren
        Reproduktionszahl ist nach Untersuchungen des Imperial College London um einen Faktor von rund 1.7 höher und sie
        kann daher zu einer beschleunigten lokalen Ausbreitung führen. Die Visualisierung basiert auf den dem Robert
        Koch-Institut gemeldeten Fällen.
      </p>
      <p>
        <b>Disclaimer:</b> This website and its contents herein, including all data and analysis are provided to the
        public strictly for educational and academic research purposes. Reliance on the Website for medical guidance or
        use of the Website in commerce is strictly prohibited.
      </p>
    </div>
  );
}
