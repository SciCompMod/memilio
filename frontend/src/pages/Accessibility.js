import React, {useEffect} from 'react';
import {Link} from 'react-router-dom';

export default function Accessibility() {
  useEffect(() => {
    document.title = document.title = `Erklärung zur Barrierefreiheit`;
  });

  return (
    <div className="accessibility">
      <Link tabIndex="1" titel="Zurück zur Hauptseite" to="/">
        Zurück zur Hauptseite
      </Link>
      <h2 className="mt-2">Erklärung zur Barrierefreiheit</h2>
      <p>
        Diese Erklärung zur Barrierefreiheit gilt für die unter hpcvscorona.dlr.de veröffentlichte Website des Deutschen
        Zentrums für Luft und Raumfahrt (DLR).
      </p>
      <p>
        Als öffentliche Stelle des Bundes sind wir bemüht, unsere Websites und mobilen Anwendungen im Einklang mit den
        Bestimmungen des Behindertengleichstellungsgesetzes des Bundes (BGG) sowie der
        Barrierefreien-Informationstechnik-Verordnung (BITV 2.0) zur Umsetzung der Richtlinie (EU) 2016/2102
        barrierefrei zugänglich zu machen.
      </p>
      <h3>Stand der Vereinbarkeit mit den Anforderungen</h3>
      <p>
        Die Anforderungen der Barrierefreiheit ergeben sich aus §§ 3 Absätze 1 bis 4 und 4 der BITV 2.0, die auf der
        Grundlage von § 12d BGG erlassen wurde.
      </p>
      <p>Die Website ist mit den zuvor genannten Anforderungen derzeit nicht vollständig vereinbar.</p>
      <p>
        Die Website wird aktuell hinsichtlich der Barrierefreiheit geprüft, daher ist es möglich, dass noch nicht alle
        Anforderungen der Barrierefreiheit vollumfänglich erfüllt sind. Wir bemühen uns, alle festgestellten Probleme
        der Zugänglichkeit zu beheben.
      </p>
      <h3>Erstellung dieser Erklärung zur Barrierefreiheit</h3>
      <p>Diese Erklärung wurde am 13. Januar 2021 erstellt.</p>
      <h3>Barrieren melden: Kontakt zu den Feedback Ansprechpartnern</h3>
      <p>
        Sie möchten uns bestehende Barrieren mitteilen oder Informationen zur Umsetzung der Barrierefreiheit erfragen?{' '}
        <br />
        Für Ihr Feedback sowie alle weiteren Informationen melden Sie sich gerne bei uns unter{' '}
        <a href="mailto:barrierefreiheit@dlr.de" title="E-Mail an: barrierefreiheit@dlr.de">
          barrierefreiheit@dlr.de
        </a>{' '}
        .
      </p>
      <h2>Schlichtungsverfahren</h2>
      <p>
        Wenn auch nach Ihrem Feedback an den oben genannten Kontakt keine zufriedenstellende Lösung gefunden wurde,
        können Sie sich an die Schlichtungsstelle nach § 16 BGG wenden. Die Schlichtungsstelle BGG hat die Aufgabe, bei
        Konflikten zum Thema Barrierefreiheit zwischen Menschen mit Behinderungen und öffentlichen Stellen des Bundes
        eine außergerichtliche Streitbeilegung zu unterstützen. Das Schlichtungsverfahren ist kostenlos. Es muss kein
        Rechtsbeistand eingeschaltet werden. Weitere Informationen zum Schlichtungsverfahren und den Möglichkeiten der
        Antragstellung erhalten Sie unter:{' '}
        <a
          href="https://www.schlichtungsstelle-bgg.de"
          target="_blank"
          rel="noopener noreferrer"
          title="Externer Link Schlichtungsstelle BGG"
        >
          www.schlichtungsstelle-bgg.de
        </a>
      </p>
      <p>
        Direkt kontaktieren können Sie die Schlichtungsstelle BGG unter:{' '}
        <a
          href="mailto:info@schlichtungsstelle-bgg.de"
          target="_blank"
          rel="noopener noreferrer"
          title="Externer Link Emailadresse Schlichtungsstelle BGG"
        >
          info@schlichtungsstelle-bgg.de
        </a>
      </p>
    </div>
  );
}
