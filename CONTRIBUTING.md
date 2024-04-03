# Contribuyendo

¡Las contribuciones son bienvenidas y se agradecen mucho! cada poquito
ayuda y siempre se dará crédito a quien colabore en este proyecto.

## Tipos de contribuciones

### Informar errores

Si está informando un error, incluya:

* Nombre y versión de su sistema operativo.
* Cualquier detalle sobre su configuración local que pueda ser útil para solucionar problemas.
* Pasos detallados para reproducir el error.

### Corrección de errores

Busque errores en *GitHub issues*. Cualquier entrada etiquetada con "bug" y/o "help wanted"
está abierto a quien quiera colaborar.

### Implementar nuevas funcionalidades

Busque entradas en los *GitHub issues* que estén etiquetadas como "enhancement"
y/o "help wanted".

### Escribir documentación

¡Nunca es posible tener suficiente documentación! Por favor siéntase libre de contribuir a cualquier
parte de la documentación, como los documentos oficiales, cadenas de documentos o incluso
en la web en publicaciones de blogs, artículos y demás.

### Enviar comentarios

Si está proponiendo una nueva característica:

* Explicar detalladamente cómo funcionaría.
* Mantenga el alcance lo más limitado posible, para que sea más fácil de implementar.
* Recuerde que este es un proyecto impulsado por voluntarios, y que las aportaciones
   son bienvenidas :)

## ¡Comenzar a contribuir!

¿Listo para contribuir? A continuación se explica cómo configurar `pytc2` para el desarrollo local.

1. Descargue una copia de `pytc2` localmente.
2. Instale `pytc2` usando `poetry`:

     ```console
     $ poetry install
     ```

3. Utilice `git` (o similar) para crear una rama para desarrollo local y realizar sus cambios:

     ```console
     $ git checkout -b name-of-your-bugfix-or-feature
     ```

4. Cuando haya terminado de realizar cambios, verifique que sus cambios cumplan con los requisitos de formato del código y pase las pruebas.

5. Haga un *commit* y abra un *pull-request*.

## Directrices para pull-request's

Antes de enviar un *pull-request*, verifique que cumpla con estas pautas:

1. El *pull-request* debe incluir pruebas adicionales, si corresponde. Revise la carpeta *test* y el formato de las pruebas incluidas.
2. Si el *pull-request* agrega funcionalidad, los documentos deben actualizarse. Revise la carpeta *docs* para ver el formato y estilo de la documentación.
3. El *pull-request* debería funcionar para todos los sistemas operativos y versiones de Python actualmente compatibles.

## Código de conducta

Tenga en cuenta que el proyecto `pytc2` se lanza con un *Código de conducta*. Al contribuir a este proyecto, usted acepta cumplir con sus términos.

