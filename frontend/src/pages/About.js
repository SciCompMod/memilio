import React from 'react';
import {Link} from 'react-router-dom';

export default function About() {
  return (
    <div className="about">
      <Link to="/">Zurück</Link>
      <h1 className="fontcolor0">HPC against Corona</h1>
      <p>
        <b>Hochaufgelöste Simulation der COVID19-Ausbreitung dank High-Performance Computing</b>
      </p>
      <p>
        Die Corona-Pandemie stellt nicht nur jeden einzelnen Menschen vor viele Fragen. Auch die Politik muss
        Entscheidungen über sinnvolle und vertretbare Gegenmaßnahmen und Eindämmungsversuche treffen. Für politische und
        wirtschaftliche Entscheidungen werden fundierte Informationen aus Bereichen der Virologie und Epidemiologie
        benötigt. Darüber hinaus kann aber auch die Expertise in den Bereichen der Mathematik und Informatik, speziell
        der numerischen Simulation einen wichtigen Beitrag leisten.
      </p>
      <p>
        Das Institut für Softwaretechnologie entwickelt in Zusammenarbeit mit dem Helmholtz-zentrum für
        Infektionsforschung ein umfassendes Softwarepaket, welches das COVID19-Infektionsgeschehen per Simulation
        darstellt. Neben einer hochaufgelösten Modellierung der Infektionsausbreitung steht im Mittelpunkt die Wirkung
        der Gegenmaßnahmen zu bewerten. Teil des Softwarepaktes ist ein Online-Tool, welches der Öffentlichkeit zu
        Verfügung gestellt wird. Das Tool soll zum einen der Politik dienen, um sich bei Entscheidungen zu weiteren
        Maßnahmen auf fundierten wissenschaftlichen Daten und Hochrechnungen zu stützen. Zum anderen soll das intuitive
        User-Interface helfen, den Verlauf der Infektionsketten auch der Bevölkerung erkennbar zu machen.
      </p>
      <p>
        <b>Relevanz von mathematischer Modellierung bei Pandemien</b>
      </p>
      <p>
        Die Komplexität der Pandemiesimulation steigt, je höher aufgelöst das Infektionsgeschehen verfolgt werden soll.
        Ein komplexes Modell ist notwendig, um verlässliche Prognosen und Auswertungen erstellen zu können. Wie in der
        kürzeren Vergangenheit geschehen, kommt es immer wieder zu lokalen Ausbrüchen der Krankheit. Diese sollen im
        besten Fall vorhergesagt und eingegrenzt werden.
      </p>
      <p>
        Mit steigender Komplexität sind konventionelle Rechner schnell überfordert, und die Rechenkraft von
        Hochleistungsrechnern wird benötigt. High-Performance Computing erlaubt, eine Vielzahl geografischer,
        demografischer und zeitlicher Faktoren im Modell zu berücksichtigen.
      </p>
      <p>
        <b>Hochaufgelöste Modellierung der Infektionsketten</b>
      </p>
      <p>
        Zusammen mit der Abteilung System Immunologie des Helmholtz-Zentrums für Infektionsforschung arbeitet das
        Institut für Softwaretechnologie an einer räumlich und demografisch hochaufgelösten Simulation der
        Corona-Pandemie für Deutschland. Auf Stadt- und Landkreisebene fließen mehrere Faktoren in die Simulation ein,
        um so unter anderem eine Bewertung der Gegenmaßnahmen auf kleinster lokaler Ebene zu ermöglichen. Die Simulation
        und die daraus resultierenden Prognosen können ebenfalls auf bundesweiter Ebene abgebildet werden.
      </p>
      <p>
        Zur Modellierung der Infektionsketten werden sowohl differentialgleichungsbasierte als auch agentenbasierte
        Modelle herangezogen. Dabei werden in der Modellierung des Krankheitsverlaufes Zustände wie „Infiziert“,
        „Trägerstatus“, „Hospitalisierung“, „Gesund“ oder „Immun“ definiert. Die möglichst altersgerechte
        Parameterschätzung in eben jene Zustände zu fallen, basiert auf großen Datenmengen verschiedener Quellen.
      </p>
      <p>
        Darüber hinaus beachtet die Simulation geografische Heterogenität. Neben der Modellierung eines
        Infektionsgeschehens auf Stadt- oder Landkreisebene fließt eine realistische Mobilität zwischen verschiedenen
        Regionen in die Simulation ein. So kann zum Beispiel der Einfluss von Pendlerverkehr untersucht werden. Die
        geografische Ausbreitung nach lokalen „Superspreading-Events“ wird ebenfalls analysiert und soll möglichst
        frühzeitig erkannt werden, um neue Infektionsketten zu verhindern.
      </p>
      <p>
        <b>Projektpartner</b>
      </p>
      <ul>
        <li>DLR-Institut für Softwaretechnologie</li>
        <li>Helmholtz-Zentrum für Infektionsforschung, Abteilung System Immunologie</li>
      </ul>
    </div>
  );
}
