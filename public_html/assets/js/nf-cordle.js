const targetPipelines = [
    ["fastqc", "cutadapt","fastqc"]
]
const FLIP_ANIMATION_DURATION = 500
const DANCE_ANIMATION_DURATION = 500
// const keyboard = document.querySelector("[data-keyboard]")
// const guessGrid = document.querySelector("[data-guess-grid]")
// const offsetFromDate = new Date(2022, 0, 1)
// const msOffset = Date.now() - offsetFromDate
// const dayOffset = msOffset / 1000 / 60 / 60 / 24
const targetPipeline = targetPipelines[0]

// startInteraction()

// function startInteraction() {
//     document.addEventListener("click", handleMouseClick)
//     document.addEventListener("keydown", handleKeyPress)
// }

// function stopInteraction() {
//     document.removeEventListener("click", handleMouseClick)
//     document.removeEventListener("keydown", handleKeyPress)
// }


// function handleKeyPress(e) {
//     if (e.key === "Enter") {
//         submitGuess()
//         return
//     }

//     if (e.key === "Backspace" || e.key === "Delete") {
//         deleteKey()
//         return
//     }

//     if (e.key.match(/^[a-z]$/)) {
//         pressKey(e.key)
//         return
//     }
// }

// function pressKey(key) {
//     const activeTiles = getActiveTiles()
//     if (activeTiles.length >= WORD_LENGTH) return
//     const nextTile = guessGrid.querySelector(":not([data-letter])")
//     nextTile.dataset.letter = key.toLowerCase()
//     nextTile.textContent = key
//     nextTile.dataset.state = "active"
// }

// function deleteKey() {
//     const activeTiles = getActiveTiles()
//     const lastTile = activeTiles[activeTiles.length - 1]
//     if (lastTile == null) return
//     lastTile.textContent = ""
//     delete lastTile.dataset.state
//     delete lastTile.dataset.letter
// }
$(document).ready(function () {
// $('.guess-tile').select2();
    $(".submit-guess").on('click', function () {
        submitGuess();
    });

});
var modules = current_pipeline_modules.split(";");
function submitGuess() {
    var guess = $('.guess-tile').map(function (i, el) { 
        return $(el).val() ? $(el) : null; });

    if (guess.length !== modules.length) {
        $(".alert-container").html('<div class="alert alert-danger">Not enough modules</div>');
        shakeTiles(guess)
        return
    }

    $('.guess-tile').map(function (i, el) { flipTile(el,i,guess)});
}

function flipTile(tile, index, array) {
    var module_guess = $(tile).val();
    setTimeout(() => {
        $(tile).addClass("flip")
    }, (index * FLIP_ANIMATION_DURATION) / 2)

    $(tile).on(
        "transitionend webkitTransitionEnd oTransitionEnd otransitionend MSTransitionEnd",
        function() {
            $(tile).removeClass("flip")
            if (modules[index] === module_guess) {
                // tile.dataset.state = "correct"
                $(tile).addClass("correct")
            } else if (modules.includes(module_guess)) {
                // tile.dataset.state = "wrong-location"
                $(tile).addClass("wrong-location")
            } else {
                // tile.dataset.state = "wrong"
                $(tile).addClass("wrong")
            }
            if (index === array.length - 1) {
                $(tile).on(
                    "transitionend webkitTransitionEnd oTransitionEnd otransitionend MSTransitionEnd",
                    () => {
                        checkWinLose(modules, array)
                    }
                )
            }
        }
    )
}

function getActiveTiles() {
    return guessGrid.querySelectorAll('.guess-tile')
}

function showAlert(message, duration = 1000) {
    const alert = document.createElement("div")
    alert.textContent = message
    $(alert).addClass("alert")
    alertContainer.prepend(alert)
    if (duration == null) return

    setTimeout(() => {
        $(alert).addClass("hide")
        alert.addEventListener("transitionend", () => {
            alert.remove()
        })
    }, duration)
}

function shakeTiles(tiles) {
    tiles.map((i, tile) => {
        $(tile).addClass("shake")
        $(tile).on(
            "animationend webkitAnimationEnd oAnimationEnd",
            () => {
                $(tile).removeClass("shake")
            }
        )
    })
}

function checkWinLose(guess, tiles) {
    if (guess === modules) {
        $(".alert-container").html('<div class="alert alert-success">Correct!</div>');
        danceTiles(tiles)
        return
    }
}

function danceTiles(tiles) {
    tiles.map((tile, index) => {
        setTimeout(() => {
            $(tile).addClass("dance")
            $(tile).on(
                "transitionend webkitTransitionEnd oTransitionEnd otransitionend MSTransitionEnd",
                function () {
                    $(tile).removeClass("dance")
                }
            )
        }, (index * DANCE_ANIMATION_DURATION) / 5)
    })
}