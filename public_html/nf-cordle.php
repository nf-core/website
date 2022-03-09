<?php
$title = 'nf-cordle';
$subtitle = 'Guess the nf-core pipeline';
$config = parse_ini_file("../config.ini");
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);

$config = parse_ini_file("../config.ini");
$conn = mysqli_connect($config['host'], $config['username'], $config['password'], $config['dbname'], $config['port']);

// get all pipelines
$sql = "SELECT * FROM nfcore_pipelines WHERE pipeline_type='pipelines' ORDER BY LOWER(name)";
$pipelines = [];
if ($result = mysqli_query($conn, $sql)) {
    if (mysqli_num_rows($result) > 0) {
        while ($row = mysqli_fetch_array($result)) {
            $pipelines[] = $row['name'];
        }
        // Free result set
        mysqli_free_result($result);
    } else {
        echo "Oops! Something went wrong. Please try again later.";
    }
}
// generate a new number every twelve hours
$date1 = date_create("2022-03-08T03:00:00");
$date2 = date_create("now");

$offset = floor(($date2->getTimestamp() - $date1->getTimestamp()) / (60 * 60 * 12));

srand($offset); // set random number seed changing every 12 hours

$target_word = $pipelines[rand(0, count($pipelines))];
$word_length = strlen($target_word);
include('../includes/header.php');
?>
<style>
    .keyboard {
        display: grid;
        grid-template-columns: repeat(20, minmax(auto, 1.25em));
        grid-auto-rows: 3em;
        gap: .25em;
        justify-content: center;
    }

    .key {
        font-size: inherit;
        grid-column: span 2;
        border: none;
        padding: 0;
        display: flex;
        justify-content: center;
        align-items: center;
        background-color: hsl(var(--hue, 200),
                var(--saturation, 1%),
                calc(var(--lightness-offset, 0%) + var(--lightness, 51%)));
        text-transform: uppercase;
        border-radius: .25em;
        cursor: pointer;
        user-select: none;
    }

    .key.large {
        grid-column: span 3;
    }

    .key>svg {
        width: 1.75em;
        height: 1.75em;
    }

    .key:hover,
    .key:focus {
        --lightness-offset: 10%;
    }

    .key.wrong {
        --lightness: 23%;
    }

    .key.wrong-location {
        --hue: 49;
        --saturation: 51%;
        --lightness: 47%;
    }

    .key.correct {
        --hue: 115;
        --saturation: 29%;
        --lightness: 43%;
    }

    .guess-grid {
        display: grid;
        justify-content: center;
        align-content: center;
        flex-grow: 1;
        grid-template-columns: repeat(<?php echo $word_length; ?>, 4em);
        grid-template-rows: repeat(calc(<?php echo $word_length + 1; ?>), 4em);
        gap: .25em;
        margin-bottom: 1em;
    }

    .tile {
        font-size: 2em;
        color: white;
        border: .05em solid hsl(240, 2%, 23%);
        text-transform: uppercase;
        font-weight: bold;
        display: flex;
        justify-content: center;
        align-items: center;
        user-select: none;
        transition: transform 250ms linear;
    }

    .tile[data-state="active"] {
        border-color: hsl(200, 1%, 34%);
    }

    .tile[data-state="wrong"] {
        border: none;
        background-color: hsl(240, 2%, 23%);
    }

    .tile[data-state="wrong-location"] {
        border: none;
        background-color: hsl(49, 51%, 47%);
    }

    .tile[data-state="correct"] {
        border: none;
        background-color: hsl(115, 29%, 43%);
    }

    .tile.shake {
        animation: shake 250ms ease-in-out;
    }

    .tile.dance {
        animation: dance 500ms ease-in-out;
    }

    .tile.flip {
        transform: rotateX(90deg);
    }

    @keyframes shake {
        10% {
            transform: translateX(-5%);
        }

        30% {
            transform: translateX(5%);
        }

        50% {
            transform: translateX(-7.5%);
        }

        70% {
            transform: translateX(7.5%);
        }

        90% {
            transform: translateX(-5%);
        }

        100% {
            transform: translateX(0);
        }
    }

    @keyframes dance {
        20% {
            transform: translateY(-50%);
        }

        40% {
            transform: translateY(5%);
        }

        60% {
            transform: translateY(-25%);
        }

        80% {
            transform: translateY(2.5%);
        }

        90% {
            transform: translateY(-5%);
        }

        100% {
            transform: translateY(0);
        }
    }

    .alert-container {
        z-index: 1;
        display: flex;
        flex-direction: column;
        align-items: center;
    }

    .alert {
        pointer-events: none;
        padding: .75em;
        border-radius: .25em;
        opacity: 1;
        transition: opacity 500ms ease-in-out;
        margin-bottom: .5em;
    }

    .alert:last-child {
        margin-bottom: 0;
    }

    .alert.hide {
        opacity: 0;
    }

    .bg-confetti-animated {
        background-repeat: repeat-x;
        background-position: top -10px center;
        background-image: url("data:image/svg+xml,%3Csvg width='600' height='90' viewBox='0 0 600 90' fill='none' xmlns='http://www.w3.org/2000/svg'%3E%3Crect x='42' y='-10' width='6' height='10'/%3E%3Crect x='84' y='-10' width='6' height='10'/%3E%3Crect x='126' y='-13' width='5' height='13'/%3E%3Crect x='168' y='-13' width='5' height='13'/%3E%3Crect x='210' y='-10' width='6' height='10'/%3E%3Crect x='252' y='-13' width='5' height='13'/%3E%3Crect x='294' y='-10' width='6' height='10'/%3E%3Crect x='336' y='-13' width='5' height='13'/%3E%3Crect x='378' y='-13' width='5' height='13'/%3E%3Crect x='420' y='-10' width='6' height='10'/%3E%3Crect x='462' y='-10' width='6' height='10'/%3E%3Crect x='504' y='-13' width='5' height='13'/%3E%3Crect x='546' y='-10' width='6' height='10'/%3E%3Cstyle type='text/css'%3E rect %7B opacity: 0; %7D rect:nth-child(1) %7B transform-origin: 45px 5px; transform: rotate(-145deg); animation: blast 700ms infinite ease-out; animation-delay: 88ms; animation-duration: 631ms; %7D rect:nth-child(2) %7B transform-origin: 87px 5px; transform: rotate(164deg); animation: blast 700ms infinite ease-out; animation-delay: 131ms; animation-duration: 442ms; %7D rect:nth-child(3) %7B transform-origin: 128px 6px; transform: rotate(4deg); animation: blast 700ms infinite ease-out; animation-delay: 92ms; animation-duration: 662ms; %7D rect:nth-child(4) %7B transform-origin: 170px 6px; transform: rotate(-175deg); animation: blast 700ms infinite ease-out; animation-delay: 17ms; animation-duration: 593ms; %7D rect:nth-child(5) %7B transform-origin: 213px 5px; transform: rotate(-97deg); animation: blast 700ms infinite ease-out; animation-delay: 122ms; animation-duration: 476ms; %7D rect:nth-child(6) %7B transform-origin: 255px 6px; transform: rotate(57deg); animation: blast 700ms infinite ease-out; animation-delay: 271ms; animation-duration: 381ms; %7D rect:nth-child(7) %7B transform-origin: 297px 5px; transform: rotate(-46deg); animation: blast 700ms infinite ease-out; animation-delay: 131ms; animation-duration: 619ms; %7D rect:nth-child(8) %7B transform-origin: 338px 6px; transform: rotate(-65deg); animation: blast 700ms infinite ease-out; animation-delay: 85ms; animation-duration: 668ms; %7D rect:nth-child(9) %7B transform-origin: 380px 6px; transform: rotate(13deg); animation: blast 700ms infinite ease-out; animation-delay: 128ms; animation-duration: 377ms; %7D rect:nth-child(10) %7B transform-origin: 423px 5px; transform: rotate(176deg); animation: blast 700ms infinite ease-out; animation-delay: 311ms; animation-duration: 508ms; %7D rect:nth-child(11) %7B transform-origin: 465px 5px; transform: rotate(108deg); animation: blast 700ms infinite ease-out; animation-delay: 108ms; animation-duration: 595ms; %7D rect:nth-child(12) %7B transform-origin: 506px 6px; transform: rotate(62deg); animation: blast 700ms infinite ease-out; animation-delay: 105ms; animation-duration: 375ms; %7D rect:nth-child(13) %7B transform-origin: 549px 5px; transform: rotate(16deg); animation: blast 700ms infinite ease-out; animation-delay: 149ms; animation-duration: 491ms; %7D rect:nth-child(odd) %7B fill: %2365BB5C; %7D rect:nth-child(even) %7B z-index: 1; fill: %2333AAFF; %7D rect:nth-child(4n) %7B animation-duration: 1400ms; fill: %23F23B14; %7D rect:nth-child(3n) %7B animation-duration: 1750ms; animation-delay: 700ms; %7D rect:nth-child(4n-7) %7B fill: %232A2F6A; %7D rect:nth-child(6n) %7B fill: %23FBBA23; %7D @keyframes blast %7B from %7B opacity: 0; %7D 20%25 %7B opacity: 1; %7D to %7B transform: translateY(90px); %7D %7D %3C/style%3E%3C/svg%3E%0A");

        @media (prefers-reduced-motion) {
            background-image: url("data:image/svg+xml,%3Csvg width='574' height='60' viewBox='0 0 574 60' fill='none' xmlns='http://www.w3.org/2000/svg'%3E%3Crect opacity='0.8' x='27.1224' y='20.0458' width='5' height='13' transform='rotate(-139 27.1224 20.0458)' fill='%23F23B14'/%3E%3Crect opacity='0.8' x='118.478' y='7.00201' width='5' height='13' transform='rotate(-38.8114 118.478 7.00201)' fill='%23FBBA23'/%3E%3Crect opacity='0.8' x='504.616' y='25.4479' width='5' height='13' transform='rotate(-60.2734 504.616 25.4479)' fill='%23F23B14'/%3E%3Crect opacity='0.6' x='538.983' y='45.555' width='5' height='13' transform='rotate(16.7826 538.983 45.555)' fill='%232A2F6A'/%3E%3Crect opacity='0.3' x='470.322' y='2.63625' width='5' height='13' transform='rotate(11.295 470.322 2.63625)' fill='%2333AAFF'/%3E%3Crect opacity='0.3' x='190.295' y='4.58138' width='5' height='13' transform='rotate(27.5954 190.295 4.58138)' fill='%23F23B14'/%3E%3Crect opacity='0.8' x='234.303' y='16.3233' width='5' height='13' transform='rotate(-41.8233 234.303 16.3233)' fill='%2365BB5C'/%3E%3Crect opacity='0.6' x='369.702' y='40.9875' width='5' height='13' transform='rotate(-56.419 369.702 40.9875)' fill='%2333AAFF'/%3E%3Crect opacity='0.3' x='402.121' y='31.0848' width='5' height='13' transform='rotate(-17.9234 402.121 31.0848)' fill='%23F23B14'/%3E%3Crect opacity='0.6' x='200.316' y='31.9328' width='5' height='13' transform='rotate(-15.8896 200.316 31.9328)' fill='%232A2F6A'/%3E%3Crect opacity='0.6' x='69.6745' y='23.4725' width='6' height='10' transform='rotate(70.0266 69.6745 23.4725)' fill='%2365BB5C'/%3E%3Crect opacity='0.6' x='291.945' y='7.16931' width='6' height='10' transform='rotate(30.4258 291.945 7.16931)' fill='%23FBBA23'/%3E%3Crect opacity='0.3' x='33.7754' y='38.2208' width='6' height='10' transform='rotate(38.6056 33.7754 38.2208)' fill='%23FBBA23'/%3E%3Crect opacity='0.8' x='109.752' y='31.1743' width='6' height='10' transform='rotate(28.5296 109.752 31.1743)' fill='%2333AAFF'/%3E%3Crect opacity='0.3' x='278.081' y='37.8695' width='6' height='10' transform='rotate(-26.5651 278.081 37.8695)' fill='%23F23B14'/%3E%3Crect opacity='0.8' x='416.294' y='11.5573' width='6' height='10' transform='rotate(-22.8498 416.294 11.5573)' fill='%23FBBA23'/%3E%3Crect opacity='0.3' x='354.667' y='9.32341' width='6' height='10' transform='rotate(17.7506 354.667 9.32341)' fill='%232A2F6A'/%3E%3Crect opacity='0.8' x='532.404' y='16.6372' width='6' height='10' transform='rotate(-75.3432 532.404 16.6372)' fill='%23FBBA23'/%3E%3Crect opacity='0.6' x='460.463' y='39.3557' width='6' height='10' transform='rotate(45.4982 460.463 39.3557)' fill='%2365BB5C'/%3E%3C/svg%3E");
        }
    }
</style>
<h1>nf-cordle/pipelines</h1>
<div class="alert-container m-auto text-center" data-alert-container></div>
<div data-guess-grid class="guess-grid ">
    <?php foreach (range(0, $word_length) as $rowindex) : ?>
        <!-- <div class=" d-flex flex-row"> -->
        <?php foreach (range(0, $word_length - 1) as $index) : ?>
            <div class="tile"></div>
        <?php endforeach; ?>
        <!-- </div> -->
    <?php endforeach; ?>
</div>
<div data-keyboard class="keyboard">
    <button class="key" data-key="Q">Q</button>
    <button class="key" data-key="W">W</button>
    <button class="key" data-key="E">E</button>
    <button class="key" data-key="R">R</button>
    <button class="key" data-key="T">T</button>
    <button class="key" data-key="Y">Y</button>
    <button class="key" data-key="U">U</button>
    <button class="key" data-key="I">I</button>
    <button class="key" data-key="O">O</button>
    <button class="key" data-key="P">P</button>
    <div class="space"></div>
    <button class="key" data-key="A">A</button>
    <button class="key" data-key="S">S</button>
    <button class="key" data-key="D">D</button>
    <button class="key" data-key="F">F</button>
    <button class="key" data-key="G">G</button>
    <button class="key" data-key="H">H</button>
    <button class="key" data-key="J">J</button>
    <button class="key" data-key="K">K</button>
    <button class="key" data-key="L">L</button>
    <div class="space"></div>
    <button data-enter class="key large">Enter</button>
    <button class="key" data-key="Z">Z</button>
    <button class="key" data-key="X">X</button>
    <button class="key" data-key="C">C</button>
    <button class="key" data-key="V">V</button>
    <button class="key" data-key="B">B</button>
    <button class="key" data-key="N">N</button>
    <button class="key" data-key="M">M</button>
    <button data-delete class="key large">
        <svg xmlns="http://www.w3.org/2000/svg" height="24" viewBox="0 0 24 24" width="24">
            <path fill="var(--color-tone-1)" d="M22 3H7c-.69 0-1.23.35-1.59.88L0 12l5.41 8.11c.36.53.9.89 1.59.89h15c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2zm0 16H7.07L2.4 12l4.66-7H22v14zm-11.59-2L14 13.41 17.59 17 19 15.59 15.41 12 19 8.41 17.59 7 14 10.59 10.41 7 9 8.41 12.59 12 9 15.59z"></path>
        </svg>
    </button>
</div>
<div class="toast congrats" role="alert" aria-live="assertive" aria-atomic="true">
    <div class="toast-header">
        <img src="/assets/img/logo/nf-core-logo-square.png" class="rounded me-2" alt="">
        <strong class="me-auto">Result copied!</strong>
        <button type="button" class="ms-2 mb-1 btn-close" data-bs-dismiss="toast" aria-label="Close"></button>
    </div>
    <div class="toast-body">
        Results copied to your clipboard. Don't be shy and humblebrag with them!
    </div>
</div>
<script type="text/javascript">
    const targetWords = <?php echo json_encode($pipelines); ?>;
    const FLIP_ANIMATION_DURATION = 500;
    const DANCE_ANIMATION_DURATION = 500;
    const offsetFromDate = new Date(2022, 0, 1);
    const msOffset = Date.now() - offsetFromDate;
    const dayOffset = msOffset / 1000 / 60 / 60 / 12;
    const targetWord = <?php echo json_encode($target_word); ?>;
    const WORD_LENGTH = targetWord.length;

    loadCookie();
    startInteraction()

    function startInteraction() {
        document.addEventListener("click", handleMouseClick)
        document.addEventListener("keydown", handleKeyPress)
    }

    function stopInteraction() {
        document.removeEventListener("click", handleMouseClick)
        document.removeEventListener("keydown", handleKeyPress)
    }

    function handleMouseClick(e) {
        if (e.target.matches("[data-key]")) {
            pressKey(e.target.dataset.key)
            return
        }

        if (e.target.matches("[data-enter]")) {
            submitGuess()
            return
        }

        if (e.target.matches("[data-delete]")) {
            deleteKey()
            return
        }
    }

    function handleKeyPress(e) {
        if (e.key === "Enter") {
            submitGuess()
            return
        }

        if (e.key === "Backspace" || e.key === "Delete") {
            deleteKey()
            return
        }

        if (e.key.match(/^[a-z]$/)) {
            pressKey(e.key)
            return
        }
    }

    function pressKey(key) {
        const activeTiles = getActiveTiles()
        if (activeTiles.length >= WORD_LENGTH) return
        const nextTile = document.querySelector("[data-guess-grid]").querySelector(":not([data-letter])")
        nextTile.dataset.letter = key.toLowerCase()
        nextTile.textContent = key
        nextTile.dataset.state = "active"
    }

    function deleteKey() {
        const activeTiles = getActiveTiles()
        const lastTile = activeTiles[activeTiles.length - 1]
        if (lastTile == null) return
        lastTile.textContent = ""
        delete lastTile.dataset.state
        delete lastTile.dataset.letter
    }

    function submitGuess() {
        const activeTiles = [...getActiveTiles()]
        const remainingTiles = document.querySelector("[data-guess-grid]").querySelectorAll(":not([data-letter])");
        const guess = activeTiles.reduce((word, tile) => {
            return word + tile.dataset.letter
        }, "")

        // if (!targetWords.includes(guess)) {
        //     showAlert("Not a name of an nf-core pipeline")
        //     shakeTiles(activeTiles)
        //     return
        // }


        for (i = 0; i < WORD_LENGTH - activeTiles.length; i++) {
            remainingTiles[i].dataset.state = "wrong"
            remainingTiles[i].dataset.letter = ""
        }

        stopInteraction()
        activeTiles.forEach((...params) => flipTile(...params, guess))

    }

    function flipTile(tile, index, array, guess) {
        const letter = tile.dataset.letter
        const key = document.querySelector("[data-keyboard]").querySelector(`[data-key="${letter}"i]`)
        setTimeout(() => {
            tile.classList.add("flip")
        }, (index * FLIP_ANIMATION_DURATION) / 2)

        tile.addEventListener(
            "transitionend",
            () => {
                tile.classList.remove("flip")
                if (targetWord[index] === letter) {
                    tile.dataset.state = "correct"
                    key.classList.add("correct")
                } else if (targetWord.includes(letter)) {
                    tile.dataset.state = "wrong-location" // bug: marks as wrong-location when the same letter, which only appears once in the correct solution, is entered twice, independent if one of the letters is already in the correct position
                    key.classList.add("wrong-location")
                } else {
                    tile.dataset.state = "wrong"
                    key.classList.add("wrong")
                }

                if (index === array.length - 1) {
                    tile.addEventListener(
                        "transitionend",
                        () => {
                            startInteraction()
                            checkWinLose(guess, array)
                            updateCookie()
                        }, {
                            once: true
                        }
                    )
                }
            }, {
                once: true
            }
        )
    }

    function getActiveTiles() {
        return document.querySelector("[data-guess-grid]").querySelectorAll('[data-state="active"]')
    }

    function showAlert(message, duration = 1000) {
        const alert = document.createElement("div")
        alert.textContent = message
        alert.classList.add("alert")
        document.querySelector("[data-alert-container]").prepend(alert)
        if (duration == null) return

        setTimeout(() => {
            alert.classList.add("hide")
            alert.addEventListener("transitionend", () => {
                alert.remove()
            })
        }, duration)
    }

    function shakeTiles(tiles) {
        tiles.forEach(tile => {
            tile.classList.add("shake")
            tile.addEventListener(
                "animationend",
                () => {
                    tile.classList.remove("shake")
                }, {
                    once: true
                }
            )
        })
    }

    function checkWinLose(guess, tiles) {
        if (guess === targetWord) {
            danceTiles(tiles)
            document.querySelector("[data-guess-grid]").classList.add("bg-confetti-animated")
            create_share_text()
            stopInteraction()
            return
        }

        const remainingTiles = document.querySelector("[data-guess-grid]").querySelectorAll(":not([data-letter])")
        if (remainingTiles.length === 0) {
            showAlert(targetWord.toUpperCase(), null)
            create_share_text()
            stopInteraction()
        }
    }

    function create_share_text() {
        var sharetext = "nf-cordle <?php echo $offset; ?>\n";
        const tiles = document.querySelectorAll("[data-guess-grid] [data-state]");
        const tilemap = Array.from(tiles).map(tile => {
            switch (tile.dataset.state) {
                case "correct":
                    return "🟩";
                case "wrong":
                    return "⬛️";
                case "wrong-location":
                    return "🟨";
                default:
                    return " ";
            }
        });
        sharetext += tilemap.join('').replace(new RegExp(`.{${WORD_LENGTH*2}}`, 'g'), '$&\n');
        navigator.clipboard.writeText(sharetext);
        $(".congrats").toast("show");
    }

    function updateCookie() {
        var cookie_state = [];
        var cookie_letter = [];
        for (let tiles of document.querySelectorAll("[data-guess-grid] [data-state]")) {
            const letter = tiles.dataset.letter
            const state = tiles.dataset.state
            if (letter == null) continue
            if (state == null) continue
            cookie_state.push(state)
            cookie_letter.push(letter)
        };
        // write to cookie
        document.cookie = "nf-cordle-state=" + JSON.stringify(cookie_state) + "; expires=Fri, 31 Dec 9999 23:59:59 GMT; path=/";
        document.cookie = "nf-cordle-letter=" + JSON.stringify(cookie_letter) + "; expires=Fri, 31 Dec 9999 23:59:59 GMT; path=/";
    }

    function loadCookie() {
        var cookie_state = [];
        var cookie_letter = [];
        if (document.cookie.indexOf("nf-cordle-state") != -1) {
            cookie_state = JSON.parse(document.cookie.split("nf-cordle-state=")[1].split(";")[0]);
        }
        if (document.cookie.indexOf("nf-cordle-letter") != -1) {
            cookie_letter = JSON.parse(document.cookie.split("nf-cordle-letter=")[1].split(";")[0]);
        }
        for (let i = 0; i < cookie_state.length; i++) {
            const tile = document.querySelectorAll("[data-guess-grid] .tile")[i]
            tile.dataset.state = cookie_state[i]
            tile.dataset.letter = cookie_letter[i]
            tile.textContent = cookie_letter[i]
        }
    }

    function danceTiles(tiles) {
        tiles.forEach((tile, index) => {
            setTimeout(() => {
                tile.classList.add("dance")
                tile.addEventListener(
                    "animationend",
                    () => {
                        tile.classList.remove("dance")
                    }, {
                        once: true
                    }
                )
            }, (index * DANCE_ANIMATION_DURATION) / 5)
        })
    }
</script>
<?php
include('../includes/footer.php');
